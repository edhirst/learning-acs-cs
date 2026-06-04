"""
Main training script for learning almost complex structures.
"""

import os

# torch.matrix_exp has no MPS implementation; enable the CPU fallback for
# Apple-Silicon users who pick fp32 (and therefore default to MPS).  Must be
# set before torch is imported.
os.environ.setdefault('PYTORCH_ENABLE_MPS_FALLBACK', '1')

import torch
import torch.optim as optim
import yaml
import json
import argparse
import sys
import functools
from pathlib import Path
from datetime import datetime
import secrets
import numpy as np
import matplotlib.pyplot as plt

print = functools.partial(print, flush=True)

from src.geometry import EvenSphere
from src.geometry.sphere_product import SphereProduct
from src.geometry.hopf_complex_structure import hopf_complex_structure
from src.models import TwoPatchAlmostComplexNet, FourPatchProductNet
from src.loss import TotalLoss
from src.loss.product_loss import ProductTotalLoss
from src.utils.chern import compute_first_chern_number


class FractionalCosineAnnealingWarmRestarts(torch.optim.lr_scheduler.LRScheduler):
    """Cosine annealing with warm restarts supporting non-integer T_mult.

    Cycle lengths are int(T_0 * T_mult^n) for n = 0, 1, 2, … so T_mult
    may be any float >= 1.  Restart epochs are precomputed up to max_epochs.
    """

    def __init__(
        self,
        optimizer: torch.optim.Optimizer,
        T_0: int,
        T_mult: float = 2.0,
        eta_min: float = 0.0,
        max_epochs: int = 100_000,
        last_epoch: int = -1,
    ):
        if T_0 < 1:
            raise ValueError(f"T_0 must be >= 1, got {T_0}")
        if T_mult < 1.0:
            raise ValueError(f"T_mult must be >= 1.0, got {T_mult}")
        self.T_0    = T_0
        self.T_mult = T_mult
        self.eta_min = eta_min

        # Precompute cumulative restart epochs and per-cycle lengths.
        self._restarts: list[int] = []  # epoch at which each cycle starts
        t = 0
        n = 0
        while t < max_epochs:
            self._restarts.append(t)
            t += max(1, int(T_0 * T_mult ** n))
            n += 1

        super().__init__(optimizer, last_epoch=last_epoch)

    def _cycle_info(self, epoch: int) -> tuple[int, int]:
        """Return (T_cur, T_i) for the given epoch."""
        # Find which cycle we're in.
        idx = 0
        for i in range(len(self._restarts) - 1):
            if self._restarts[i] <= epoch < self._restarts[i + 1]:
                idx = i
                break
        else:
            idx = len(self._restarts) - 1
        T_cur = epoch - self._restarts[idx]
        if idx + 1 < len(self._restarts):
            T_i = self._restarts[idx + 1] - self._restarts[idx]
        else:
            T_i = max(1, int(self.T_0 * self.T_mult ** idx))
        return T_cur, T_i

    def get_lr(self) -> list[float]:  # type: ignore[override]
        T_cur, T_i = self._cycle_info(self.last_epoch)
        return [
            self.eta_min + (base_lr - self.eta_min) * (1 + np.cos(np.pi * T_cur / T_i)) / 2
            for base_lr in self.base_lrs
        ]

class Trainer:
    """Training manager for almost complex structure learning."""
    
    def __init__(self, config_path: str):
        """Initialize trainer with configuration."""
        # Load configuration
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        # Set random seed
        torch.manual_seed(self.config['seed'])
        np.random.seed(self.config['seed'])
        
        # Setup device
        if torch.cuda.is_available():
            self.device = torch.device('cuda')
        elif torch.backends.mps.is_available():
            self.device = torch.device('mps')
        else:
            self.device = torch.device('cpu')
        print(f"Using device: {self.device}")

        # Setup floating-point precision.  The canonical config field is
        # `training.precision` (integer 32 or 64).  The legacy field
        # `training.dtype` ("float32"/"float64") is honoured as a fallback
        # so older configs continue to load.  Default: 32.
        train_cfg = self.config['training']
        if 'precision' in train_cfg:
            precision = int(train_cfg['precision'])
        elif 'dtype' in train_cfg:
            precision = 64 if str(train_cfg['dtype']).endswith('64') else 32
        else:
            precision = 32
        if precision not in (32, 64):
            raise ValueError(
                f"`training.precision` must be 32 or 64, got {precision!r}"
            )
        self.dtype = torch.float64 if precision == 64 else torch.float32
        if self.dtype == torch.float64 and self.device.type == 'mps':
            print("Warning: MPS does not support float64; falling back to CPU.")
            self.device = torch.device('cpu')
        # `torch.matrix_exp` has no MPS kernel.  Running it on MPS forces a
        # CPU round-trip on every call, which dominates the training step and
        # erases (or even reverses) the fp32 matmul speedup.  Fall back to
        # CPU unless the user opts into the Taylor-series matrix_exp
        # (`network.matrix_exp_terms` set to a positive integer), which is
        # pure matmul and runs natively on MPS.
        matrix_exp_terms = self.config.get('network', {}).get('matrix_exp_terms', None)
        if self.device.type == 'mps' and (matrix_exp_terms is None or int(matrix_exp_terms) <= 0):
            print("Note: torch.matrix_exp has no MPS kernel; falling back to "
                  "CPU to avoid per-step MPS↔CPU round-trips.")
            self.device = torch.device('cpu')
        print(f"Using precision: float{precision}")
        
        # Setup geometry and model
        geo = self.config['geometry']
        self.manifold_type = geo.get('manifold_type', 'sphere')

        if self.manifold_type == 'sphere_product':
            dim_alpha = geo['real_dimension_alpha']
            dim_beta  = geo['real_dimension_beta']
            if (dim_alpha + dim_beta) % 2 != 0:
                print(
                    f"ERROR: real_dimension_alpha + real_dimension_beta must be even for an almost "
                    f"complex structure to exist on the product; "
                    f"got {dim_alpha} + {dim_beta} = {dim_alpha + dim_beta} (odd)."
                )
                sys.exit(1)
            self.sphere_product = SphereProduct(dim_alpha, dim_beta, device=str(self.device))
            d = self.sphere_product.total_dim
            n = self.sphere_product.complex_dim
            print(f"Manifold: S^{dim_alpha} × S^{dim_beta}  (real dim {d}, complex dim {n})")

            self.model = FourPatchProductNet(
                dim_alpha=dim_alpha,
                dim_beta=dim_beta,
                hidden_layers=self.config['network']['hidden_layers'],
                activation=self.config['network']['activation'],
                use_fourier_features=self.config['network']['use_fourier_features'],
                num_fourier_features=self.config['network']['num_fourier_features'],
                fourier_scale=self.config['network']['fourier_scale'],
                initialization=self.config['network']['initialization'],
                matrix_exp_terms=self.config['network'].get('matrix_exp_terms'),
                use_residual_blocks=self.config['network'].get('use_residual_blocks', False),
            ).to(self.device).to(self.dtype)

            self.loss_fn = ProductTotalLoss(
                dim_alpha=dim_alpha,
                dim_beta=dim_beta,
                weights=self.config['loss_weights'],
                compute_nijenhuis_rank=self.config['logging'].get('compute_nijenhuis_rank', False),
            ).to(self.device).to(self.dtype)

        else:
            real_dim = geo['real_dimension']
            if not isinstance(real_dim, int) or real_dim <= 0 or real_dim % 2 != 0:
                print(
                    f"ERROR: real_dimension must be a positive even integer for manifold_type 'sphere'; "
                    f"got {real_dim!r}."
                )
                sys.exit(1)
            self.sphere = EvenSphere(
                dimension=real_dim,
                device=str(self.device),
            )
            print(f"Manifold: S^{self.sphere.dim} (complex dimension {self.sphere.complex_dim})")

            self.model = TwoPatchAlmostComplexNet(
                manifold_dim=self.sphere.dim,
                hidden_layers=self.config['network']['hidden_layers'],
                activation=self.config['network']['activation'],
                use_fourier_features=self.config['network']['use_fourier_features'],
                num_fourier_features=self.config['network']['num_fourier_features'],
                fourier_scale=self.config['network']['fourier_scale'],
                initialization=self.config['network']['initialization'],
                matrix_exp_terms=self.config['network'].get('matrix_exp_terms'),
                use_residual_blocks=self.config['network'].get('use_residual_blocks', False),
            ).to(self.device).to(self.dtype)

            self.loss_fn = TotalLoss(
                manifold_dim=self.sphere.dim,
                weights=self.config['loss_weights'],
                compute_nijenhuis_rank=self.config['logging'].get('compute_nijenhuis_rank', False),
            ).to(self.device).to(self.dtype)

        print(f"Model parameters: {sum(p.numel() for p in self.model.parameters())}")

        if (
            self.config['training'].get('equatorial_sampling', False)
            and float(self.config['loss_weights'].get('nijenhuis', 0.0)) != 0.0
        ):
            print(
                "ERROR: equatorial_sampling=true is incompatible with a non-zero nijenhuis "
                "weight. Equatorial sampling restricts points to the overlap band, so the "
                "Nijenhuis loss never sees the near-pole region where ∂J/∂x is unconstrained."
            )
            sys.exit(1)

        # Setup optimizer
        self.optimizer = self._create_optimizer()
        self.scheduler = self._create_scheduler()
        
        # Setup logging
        self.setup_logging()
        
        # Training state
        self.epoch = 0
        self.best_loss = float('inf')
        self.patience_counter = 0
        
    def _create_optimizer(self):
        """Create optimizer from config."""
        opt_name = self.config['training']['optimizer'].lower()
        lr = self.config['training']['learning_rate']
        weight_decay = self.config['training']['weight_decay']
        
        if opt_name == 'adam':
            return optim.Adam(self.model.parameters(), lr=lr, weight_decay=weight_decay)
        elif opt_name == 'adamw':
            return optim.AdamW(self.model.parameters(), lr=lr, weight_decay=weight_decay)
        elif opt_name == 'sgd':
            return optim.SGD(self.model.parameters(), lr=lr, weight_decay=weight_decay, 
                           momentum=0.9)
        else:
            raise ValueError(f"Unknown optimizer: {opt_name}")

    def _create_lbfgs_optimizer(self):
        """Create the optional post-Adam L-BFGS optimizer from config."""
        training_cfg = self.config['training']
        line_search_fn = training_cfg.get('lbfgs_line_search_fn', 'strong_wolfe')
        if isinstance(line_search_fn, str) and line_search_fn.lower() in ('none', ''):
            line_search_fn = None

        kwargs = {
            'lr': float(training_cfg.get('lbfgs_learning_rate', 1.0)),
            'max_iter': int(training_cfg.get('lbfgs_max_iter', 20)),
            'history_size': int(training_cfg.get('lbfgs_history_size', 100)),
            'line_search_fn': line_search_fn,
        }
        if 'lbfgs_max_eval' in training_cfg:
            kwargs['max_eval'] = int(training_cfg['lbfgs_max_eval'])
        if 'lbfgs_tolerance_grad' in training_cfg:
            kwargs['tolerance_grad'] = float(training_cfg['lbfgs_tolerance_grad'])
        if 'lbfgs_tolerance_change' in training_cfg:
            kwargs['tolerance_change'] = float(training_cfg['lbfgs_tolerance_change'])

        return optim.LBFGS(self.model.parameters(), **kwargs)
    
    def _create_scheduler(self):
        """Create learning rate scheduler from config."""
        sched_name = self.config['training']['scheduler'].lower()
        
        if sched_name == 'none':
            return None
        elif sched_name == 'cosine':
            params = self.config['training']['scheduler_params']
            return optim.lr_scheduler.CosineAnnealingLR(
                self.optimizer,
                T_max=params['T_max'],
                eta_min=params['eta_min']
            )
        elif sched_name == 'cosine_restarts':
            params = self.config['training']['scheduler_params']
            return FractionalCosineAnnealingWarmRestarts(
                self.optimizer,
                T_0=int(params['T_0']),
                T_mult=float(params.get('T_mult', 2)),
                eta_min=float(params['eta_min']),
                max_epochs=self.config['training']['epochs'],
            )
        elif sched_name == 'step':
            return optim.lr_scheduler.StepLR(
                self.optimizer,
                step_size=1000,
                gamma=0.5
            )
        elif sched_name == 'exponential':
            return optim.lr_scheduler.ExponentialLR(
                self.optimizer,
                gamma=0.99
            )
        else:
            raise ValueError(f"Unknown scheduler: {sched_name}")
    
    def setup_logging(self):
        """Setup output directories."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + secrets.token_hex(3)
        self.output_dir = Path(self.config['logging']['save_dir']) / timestamp
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.checkpoint_dir = self.output_dir / 'checkpoints'
        self.checkpoint_dir.mkdir(exist_ok=True)

        self.loss_data: dict[str, list | dict] = {'train': [], 'val': [], 'test': {}}

        # Save config
        with open(self.output_dir / 'config.yaml', 'w') as f:
            yaml.dump(self.config, f)

        print(f"Outputs will be saved to: {self.output_dir}")
    
    def sample_batch(self, batch_size: int):
        """Sample a batch of coordinates.

        Returns:
            For 'sphere': a single Tensor of patch-1 ball coords.
            For 'sphere_product': a tuple (x_alpha, x_beta) of ball coords.
        """
        equatorial = self.config['training'].get('equatorial_sampling', False)
        density_power = float(self.config['training'].get('density_power', 1.0))
        if self.manifold_type == 'sphere_product':
            lw = self.config['loss_weights']
            osa = float(lw.get('overlap_size_alpha', 0.1)) if equatorial else None
            osb = float(lw.get('overlap_size_beta',  0.1)) if equatorial else None
            xa, xb = self.sphere_product.sample_points(batch_size, osa, osb, density_power)
            return xa.to(device=self.device, dtype=self.dtype), xb.to(device=self.device, dtype=self.dtype)
        os_ = float(self.config['loss_weights'].get('overlap_size', 0.1)) if equatorial else None
        return self.sphere.sample_points_patch1(batch_size, density_power=density_power, overlap_size=os_).to(device=self.device, dtype=self.dtype)
    
    def _compute_train_loss_and_info(self, x) -> tuple[torch.Tensor, dict]:
        """Compute the current training loss and logging info for one batch."""
        if self.manifold_type == 'sphere_product':
            x_alpha, x_beta = x
            J_dict, x_alpha, x_beta, x_a2, x_b2 = self.model(x_alpha, x_beta)
            loss, info = self.loss_fn(J_dict, x_alpha, x_beta, x_a2, x_b2, self.model)
            if self.config['logging'].get('compute_block_diagonality', False):
                with torch.no_grad():
                    info['block_diagonality'] = self._block_diagonality(
                        J_dict, self.config['geometry']['real_dimension_alpha']
                    )
        else:
            J1, J2, x1, x2 = self.model(x)
            loss, info = self.loss_fn(J1, J2, x1, x2, self.model)

        return loss, info

    def _clip_gradients(self):
        """Apply configured gradient clipping to model parameters."""
        if self.config['advanced']['gradient_clipping'] > 0:
            torch.nn.utils.clip_grad_norm_(
                self.model.parameters(),
                self.config['advanced']['gradient_clipping']
            )

    def train_step(self, x) -> dict:
        """Perform one first-order optimizer training step."""
        self.model.train()
        self.optimizer.zero_grad()

        loss, info = self._compute_train_loss_and_info(x)

        # Backward pass
        loss.backward()
        self._clip_gradients()
        
        self.optimizer.step()
        
        return info

    def train_step_lbfgs(self, x) -> dict:
        """Perform one L-BFGS optimizer step using a fixed-batch closure."""
        self.model.train()
        info_holder: dict[str, dict] = {}

        def closure():
            self.optimizer.zero_grad()
            loss, info = self._compute_train_loss_and_info(x)
            loss.backward()
            self._clip_gradients()
            info_holder['info'] = info
            return loss

        self.optimizer.step(closure)

        # L-BFGS may evaluate the closure multiple times. Return the metrics
        # from the last evaluation, matching the loss seen by the optimizer.
        return info_holder.get('info', {})
    
    def validate(self) -> dict:
        """Run validation."""
        self.model.eval()

        num_val_points = self.config['training']['batch_size_val']
        x = self.sample_batch(num_val_points)

        with torch.no_grad():
            if self.manifold_type == 'sphere_product':
                x_alpha, x_beta = x
                J_dict, x_alpha, x_beta, x_a2, x_b2 = self.model(x_alpha, x_beta)
                _, info = self.loss_fn(J_dict, x_alpha, x_beta, x_a2, x_b2, self.model,
                                       hard_cutoff=True)
                if self.config['logging'].get('compute_block_diagonality', False):
                    info['block_diagonality'] = self._block_diagonality(
                        J_dict, self.config['geometry']['real_dimension_alpha']
                    )
            else:
                J1, J2, x1, x2 = self.model(x)
                _, info = self.loss_fn(J1, J2, x1, x2, self.model, hard_cutoff=True)

        return info
    
    def pretrain(self):
        """Warm-start each PatchSubNet by driving A(x) → 0.

        With A = 0, patch1 outputs J₀ and patch2 outputs −J₀, which is the
        correct initialisation given that the north→south chart transition is
        complex-anti-linear (J → −J).  The loss is mean(A²).
        Runs only when config['pretraining']['identity_pretraining'] is True.
        """
        pt_cfg = self.config.get('pretraining', {})
        if not pt_cfg.get('identity_pretraining', False):
            return

        pretrain_epochs = pt_cfg['pretrain_epochs']
        pretrain_points = pt_cfg['pretrain_points']

        if self.manifold_type == 'sphere_product':
            subnets = [(f'patch_{ia}{ib}', sub)
                       for (ia, ib), sub in zip(
                           [(0,0),(0,1),(1,0),(1,1)], self.model.subnets
                       )]
            print(f"\nPretraining: A→0 warm-start ({pretrain_epochs} epochs, "
                  f"{pretrain_points} points/step, 4 patches)")
        else:
            subnets = [('patch1', self.model.patch1), ('patch2', self.model.patch2)]
            print(f"\nPretraining: A→0 warm-start ({pretrain_epochs} epochs, "
                  f"{pretrain_points} points/step)")

        for subnet_name, subnet in subnets:
            opt = torch.optim.Adam(
                subnet.parameters(),
                lr=self.config['training']['learning_rate']
            )
            for epoch in range(pretrain_epochs):
                subnet.train()
                opt.zero_grad()
                if self.manifold_type == 'sphere_product':
                    # Each subnet receives the full combined coords matching
                    # the patch's α/β index pair.
                    xa, xb = self.sample_batch(pretrain_points)
                    ia, ib = [(0,0),(0,1),(1,0),(1,1)][subnets.index((subnet_name, subnet))]
                    from src.geometry.sphere import EvenSphere as _ES
                    xa_p = xa if ia == 0 else _ES.patch_change_coords(xa)
                    xb_p = xb if ib == 0 else _ES.patch_change_coords(xb)
                    x_pt = torch.cat([xa_p, xb_p], dim=1)
                else:
                    x_pt = self.sample_batch(pretrain_points)
                A = subnet.forward_raw(x_pt)              # (N, d, d)
                # Target: A = 0 so that J = J₀ (patch 1) or J = -J₀ (patch 2),
                # the correct initialisation for the two chart orientations.
                loss = torch.mean(A ** 2)
                loss.backward()
                if self.config['advanced']['gradient_clipping'] > 0:
                    torch.nn.utils.clip_grad_norm_(
                        subnet.parameters(),
                        self.config['advanced']['gradient_clipping']
                    )
                opt.step()
                if epoch % max(1, pretrain_epochs // 10) == 0:
                    print(f"  [{subnet_name}] pretrain epoch {epoch:5d} | "
                          f"mean(A²): {loss.item():.4e}")

        print("Pretraining complete.\n")

    # ------------------------------------------------------------------
    # Overlap pretraining helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _pushforward_target_sphere(
        x: torch.Tensor,
        patch_idx: int,
        J0_ref: torch.Tensor,
    ) -> torch.Tensor:
        """Target J on a single sphere that satisfies the overlap condition.

        For patch_idx = 0 (reference patch):  J*(x) = J0_ref                (constant).
        For patch_idx = 1 (non-reference):    J*(x) = F(x)^{-1} J0_ref F(x),
                                              with F = EvenSphere.patch_change_jacobian(x).

        Using the self-inverse identity F(x1) F(x2) = I (x2 = φ(x1)), this is the
        pushforward of the reference J0 expressed in the non-reference chart, so the
        overlap condition J_2 = F(x1) J_1 F(x1)^{-1} is satisfied pointwise.

        Args:
            x:         Ball coordinates in this patch's chart, shape (batch, d).
            patch_idx: 0 for reference patch, 1 for non-reference.
            J0_ref:    Canonical J0 of the reference patch, shape (d, d).

        Returns:
            Target J tensor of shape (batch, d, d).
        """
        batch = x.shape[0]
        if patch_idx == 0:
            return J0_ref.unsqueeze(0).expand(batch, -1, -1).contiguous()
        F     = EvenSphere.patch_change_jacobian(x)          # (batch, d, d)
        F_inv = torch.linalg.inv(F)                           # (batch, d, d)
        return F_inv @ J0_ref.unsqueeze(0) @ F                # (batch, d, d)

    @staticmethod
    def _pushforward_target_product(
        xa: torch.Tensor,
        xb: torch.Tensor,
        ia: int,
        ib: int,
        J0_ref: torch.Tensor,
        dim_a: int,
        dim_b: int,
    ) -> torch.Tensor:
        """Target J on a sphere product that satisfies overlap for all patch pairs.

        J*(x) = F^{-1} J0_ref F  with  F = block-diag(F_a, F_b), where
            F_a = I                                 if ia = 0, else EvenSphere.patch_change_jacobian(xa)
            F_b = I                                 if ib = 0, else EvenSphere.patch_change_jacobian(xb)

        Because F is block-diagonal, this target reduces to the factor-wise pushforward.
        The derivation verifies that it satisfies the overlap condition simultaneously
        for all six patch-pair transitions.

        Args:
            xa:     α-factor coords in chart ia, shape (batch, dim_a).
            xb:     β-factor coords in chart ib, shape (batch, dim_b).
            ia, ib: Chart indices (0 or 1) for α- and β-factors.
            J0_ref: Reference J0 at patch (0,0), shape (d, d) with d = dim_a + dim_b.
            dim_a:  Real α-factor dimension.
            dim_b:  Real β-factor dimension.

        Returns:
            Target J tensor of shape (batch, d, d).
        """
        batch = xa.shape[0]
        d     = dim_a + dim_b

        F_full = torch.zeros(batch, d, d, dtype=xa.dtype, device=xa.device)
        if ia == 0:
            F_full[:, :dim_a, :dim_a] = torch.eye(dim_a, dtype=xa.dtype, device=xa.device)
        else:
            F_full[:, :dim_a, :dim_a] = EvenSphere.patch_change_jacobian(xa)
        if ib == 0:
            F_full[:, dim_a:, dim_a:] = torch.eye(dim_b, dtype=xa.dtype, device=xa.device)
        else:
            F_full[:, dim_a:, dim_a:] = EvenSphere.patch_change_jacobian(xb)

        F_inv = torch.linalg.inv(F_full)
        return F_inv @ J0_ref.unsqueeze(0) @ F_full

    def pretrain_overlap(self):
        """Warm-start so the overlap loss starts at ~zero on the overlap band.

        Each subnet's J is supervised toward the pushforward target

            J*(x) = F^{-1} · J0_ref · F,

        where F is the Jacobian of the chart change from this patch back to the
        reference patch (patch 1 for sphere, patch (0,0) for sphere_product), and
        J0_ref is the reference patch's canonical J0.  By construction this target
        satisfies the overlap condition  J_j = F^{(j←i)} J_i (F^{(j←i)})^{-1}
        pointwise, so the overlap loss is ~zero at the end of this pretraining.

        Sampling is restricted to the equatorial annulus of each factor
        (the overlap band where F is bounded and the overlap loss is actually
        evaluated); outside this band the pushforward target is unconstrained
        and relies on network smoothness to extrapolate.

        Loss: direct Frobenius MSE between J_net(x) and J*(x).  No analytical
        inversion of the A → J exp-map is required.

        Runs only when config['pretraining']['overlap_pretraining'] is True.
        Composes cleanly with identity pretraining: if both are enabled,
        identity runs first (giving a globally smooth J ≈ ±J0 start), then
        overlap pretraining refines A away from 0 within the band.
        """
        pt_cfg = self.config.get('pretraining', {})
        if not pt_cfg.get('overlap_pretraining', False):
            return

        pretrain_epochs = pt_cfg['overlap_pretrain_epochs']
        pretrain_points = pt_cfg['overlap_pretrain_points']
        pt_lr = float(pt_cfg.get('overlap_pretrain_lr',
                                 self.config['training']['learning_rate']))
        lw = self.config['loss_weights']

        if self.manifold_type == 'sphere_product':
            os_a = float(lw.get('overlap_size_alpha', 0.1))
            os_b = float(lw.get('overlap_size_beta',  0.1))
            dim_a = self.config['geometry']['real_dimension_alpha']
            dim_b = self.config['geometry']['real_dimension_beta']

            # Reference J0 (patch (0,0)); block-diag for even-even factors,
            # full canonical J0 for odd-odd factors.  Either way, we use it
            # as-is — the formula J* = F^{-1} J0_ref F is correct in both cases.
            J0_ref = self.model.patches['patch_00'].J0.detach().to(
                device=self.device, dtype=self.dtype,
            )

            subnet_list = [
                ((ia, ib), self.model.patches[f'patch_{ia}{ib}'])
                for ia, ib in [(0, 0), (0, 1), (1, 0), (1, 1)]
            ]

            print(f"\nOverlap pretraining: J → F⁻¹ J₀ F pushforward "
                  f"({pretrain_epochs} epochs, {pretrain_points} points/step, 4 patches)")

            for (ia, ib), subnet in subnet_list:
                opt = torch.optim.Adam(subnet.parameters(), lr=pt_lr)
                for epoch in range(pretrain_epochs):
                    subnet.train()
                    opt.zero_grad()

                    # Sample from both-factor equatorial annulus.
                    xa, xb = self.sphere_product.sample_points(
                        pretrain_points, overlap_size_alpha=os_a, overlap_size_beta=os_b,
                    )
                    xa = xa.to(device=self.device, dtype=self.dtype)
                    xb = xb.to(device=self.device, dtype=self.dtype)

                    # Transform to this patch's coords.
                    xa_p = xa if ia == 0 else EvenSphere.patch_change_coords(xa)
                    xb_p = xb if ib == 0 else EvenSphere.patch_change_coords(xb)

                    # Analytical pushforward target.
                    with torch.no_grad():
                        J_target = self._pushforward_target_product(
                            xa_p, xb_p, ia, ib, J0_ref, dim_a, dim_b,
                        )

                    # Forward through this subnet.
                    x_pt  = torch.cat([xa_p, xb_p], dim=1)
                    J_net = subnet(x_pt)

                    # Mean squared Frobenius error per sample.
                    loss = ((J_net - J_target) ** 2).sum(dim=(1, 2)).mean()
                    loss.backward()
                    if self.config['advanced']['gradient_clipping'] > 0:
                        torch.nn.utils.clip_grad_norm_(
                            subnet.parameters(),
                            self.config['advanced']['gradient_clipping'],
                        )
                    opt.step()

                    if epoch % max(1, pretrain_epochs // 10) == 0:
                        print(f"  [patch_{ia}{ib}] overlap pretrain epoch {epoch:5d} | "
                              f"‖J−J*‖²_F: {loss.item():.4e}")

        else:  # sphere
            os_ = float(lw.get('overlap_size', 0.1))

            # Reference J0 = patch 1's J0 (canonical, j0_sign=+1).
            J0_ref = self.model.patch1.J0.detach().to(
                device=self.device, dtype=self.dtype,
            )

            subnet_list = [(0, self.model.patch1), (1, self.model.patch2)]

            print(f"\nOverlap pretraining: J → F⁻¹ J₀ F pushforward "
                  f"({pretrain_epochs} epochs, {pretrain_points} points/step, 2 patches)")

            for patch_idx, subnet in subnet_list:
                opt = torch.optim.Adam(subnet.parameters(), lr=pt_lr)
                for epoch in range(pretrain_epochs):
                    subnet.train()
                    opt.zero_grad()

                    # Sample from patch-1 equatorial annulus (self-symmetric under φ).
                    x1 = self.sphere.sample_points_patch1(
                        pretrain_points, overlap_size=os_,
                    ).to(device=self.device, dtype=self.dtype)

                    # Transform to this patch's coords.
                    x_p = x1 if patch_idx == 0 else EvenSphere.patch_change_coords(x1)

                    # Analytical pushforward target.
                    with torch.no_grad():
                        J_target = self._pushforward_target_sphere(x_p, patch_idx, J0_ref)

                    # Forward through this subnet.
                    J_net = subnet(x_p)

                    # Mean squared Frobenius error per sample.
                    loss = ((J_net - J_target) ** 2).sum(dim=(1, 2)).mean()
                    loss.backward()
                    if self.config['advanced']['gradient_clipping'] > 0:
                        torch.nn.utils.clip_grad_norm_(
                            subnet.parameters(),
                            self.config['advanced']['gradient_clipping'],
                        )
                    opt.step()

                    if epoch % max(1, pretrain_epochs // 10) == 0:
                        patch_name = f"patch{patch_idx + 1}"
                        print(f"  [{patch_name}] overlap pretrain epoch {epoch:5d} | "
                              f"‖J−J*‖²_F: {loss.item():.4e}")

        print("Overlap pretraining complete.\n")

    def pretrain_hopf(self):
        """Warm-start each patch subnet toward the analytic Hopf complex structure.

        The Hopf manifold construction gives an integrable complex structure on
        S^{dim_alpha} × S¹ (dim_alpha odd).  Each patch subnet is supervised
        toward J_hopf = (DF)⁻¹ J₀ DF (computed analytically in ball coords),
        driving both the Nijenhuis loss and the overlap loss toward zero.

        Loss: mean Frobenius ‖J_net(x) − J_hopf(x)‖²_F per sample.
        Sampling covers the full ball (not restricted to the equatorial band).
        Samples with ‖x_beta‖ ≥ 0.9 are filtered out to avoid near-pole
        singularities in the S¹ factor.

        Runs only when config['pretraining']['hopf_pretraining'] is True.
        Only valid for sphere_product with dim_beta = 1.
        """
        pt_cfg = self.config.get('pretraining', {})
        if not pt_cfg.get('hopf_pretraining', False):
            return

        if self.manifold_type != 'sphere_product':
            print(
                "ERROR: hopf_pretraining requires manifold_type = 'sphere_product'; "
                f"got '{self.manifold_type}'."
            )
            sys.exit(1)

        dim_beta = self.config['geometry']['real_dimension_beta']
        if dim_beta != 1:
            print(
                f"ERROR: hopf_pretraining requires real_dimension_beta = 1 "
                f"(S¹ factor); got {dim_beta}."
            )
            sys.exit(1)

        pretrain_epochs = pt_cfg['pretrain_epochs']
        pretrain_points = pt_cfg['pretrain_points']
        pt_lr = float(pt_cfg.get('hopf_pretrain_lr',
                                  self.config['training']['learning_rate']))
        dim_alpha = self.config['geometry']['real_dimension_alpha']

        subnet_list = [
            ((ia, ib), self.model.patches[f'patch_{ia}{ib}'])
            for ia, ib in [(0, 0), (0, 1), (1, 0), (1, 1)]
        ]

        print(f"\nHopf pretraining: J → (DF)⁻¹ J₀ DF  "
              f"({pretrain_epochs} epochs, {pretrain_points} points/step, 4 patches)")

        for (ia, ib), subnet in subnet_list:
            opt = torch.optim.Adam(subnet.parameters(), lr=pt_lr)
            for epoch in range(pretrain_epochs):
                subnet.train()
                opt.zero_grad()

                # Sample full-domain points.
                xa, xb = self.sphere_product.sample_points(pretrain_points)
                xa = xa.to(device=self.device, dtype=self.dtype)
                xb = xb.to(device=self.device, dtype=self.dtype)

                # Filter out samples near sphere poles where DF is ill-conditioned.
                # For pole=1 charts, patch_change_coords maps |x|≈0 → |φ(x)|≈1 (singular).
                xa2 = EvenSphere.patch_change_coords(xa)
                xb2 = EvenSphere.patch_change_coords(xb)
                mask = (xa.norm(dim=1) < 0.85) & (xb.norm(dim=1) < 0.85) \
                     & (xa2.norm(dim=1) < 0.85) & (xb2.norm(dim=1) < 0.85)
                if mask.sum() < 2:
                    continue
                xa, xb = xa[mask], xb[mask]

                # Transform to this patch's coordinates.
                xa_p = xa if ia == 0 else EvenSphere.patch_change_coords(xa)
                xb_p = xb if ib == 0 else EvenSphere.patch_change_coords(xb)

                # Analytic Hopf target for this patch.
                with torch.no_grad():
                    J_target = hopf_complex_structure(xa_p, xb_p, ia, ib)

                # Forward through this subnet.
                x_pt  = torch.cat([xa_p, xb_p], dim=1)
                J_net = subnet(x_pt)

                # Mean squared Frobenius error per sample.
                loss = ((J_net - J_target) ** 2).sum(dim=(1, 2)).mean()
                loss.backward()
                if self.config['advanced']['gradient_clipping'] > 0:
                    torch.nn.utils.clip_grad_norm_(
                        subnet.parameters(),
                        self.config['advanced']['gradient_clipping'],
                    )
                opt.step()

                if epoch % max(1, pretrain_epochs // 10) == 0:
                    print(f"  [patch_{ia}{ib}] hopf pretrain epoch {epoch:5d} | "
                          f"‖J−J_hopf‖²_F: {loss.item():.4e}")

        print("Hopf pretraining complete.\n")

    def train(self):
        """Main training loop."""
        # Identity and overlap pretraining are mutually exclusive alternatives:
        # identity drives A → 0 (so J ≡ ±J_0, Nijenhuis loss = 0 at start),
        # overlap drives J → F⁻¹J_0F  (so overlap loss ≈ 0 at start).
        # They pull A in different directions; composing them silently would
        # undo the first or give unintended intermediate states.
        pt_cfg = self.config.get('pretraining', {})
        active = [
            name for name in ('identity_pretraining', 'overlap_pretraining', 'hopf_pretraining')
            if pt_cfg.get(name, False)
        ]
        if len(active) > 1:
            print(
                f"ERROR: At most one pretraining mode may be enabled at a time; "
                f"got {active}. "
                "identity_pretraining drives A → 0 (J ≡ ±J_0), "
                "overlap_pretraining drives J → F⁻¹ J_0 F, "
                "hopf_pretraining drives J → (DF)⁻¹ J_0 DF (Hopf CS). "
                "Enable exactly one."
            )
            sys.exit(1)

        self.pretrain()
        self.pretrain_overlap()
        self.pretrain_hopf()

        print("\nStarting training...")
        print(f"Epochs: {self.config['training']['epochs']}")
        print(f"Batch size: {self.config['training']['batch_size']}")
        
        num_epochs = self.config['training']['epochs']
        use_lbfgs = self.config['training'].get('use_lbfgs_after_adam', False)
        lbfgs_epochs = int(self.config['training'].get('lbfgs_epochs', 0)) if use_lbfgs else 0
        total_epochs = num_epochs + lbfgs_epochs
        batch_size = self.config['training']['batch_size']
        log_interval = self.config['logging']['log_interval']
        checkpoint_interval = self.config['logging']['checkpoint_interval']
        if lbfgs_epochs > 0:
            print(f"L-BFGS fine-tuning epochs: {lbfgs_epochs}")
        if self.manifold_type == 'sphere_product':
            history_keys = ['total_loss', 'nij_00', 'nij_01', 'nij_10', 'nij_11', 'overlap_norm']
            if self.config['logging'].get('compute_block_diagonality', False):
                history_keys.append('block_diagonality')
        else:
            history_keys = ['total_loss', 'nij1', 'nij2', 'overlap_norm']
        history: dict[str, list] = {'epoch': [], **{k: [] for k in history_keys}}
        
        for epoch in range(total_epochs):
            using_lbfgs = epoch >= num_epochs
            if using_lbfgs and epoch == num_epochs:
                print("\nSwitching optimizer to L-BFGS...")
                self.optimizer = self._create_lbfgs_optimizer()
                self.scheduler = None

            self.epoch = epoch
            
            # Sample batch
            x = self.sample_batch(batch_size)
            
            # Training step
            train_info = self.train_step_lbfgs(x) if using_lbfgs else self.train_step(x)
            
            # Logging
            is_last_adam_epoch = lbfgs_epochs > 0 and epoch == num_epochs - 1
            is_last_epoch = (epoch == total_epochs - 1) or is_last_adam_epoch
            if epoch % log_interval == 0 or is_last_epoch:
                self._log_metrics(train_info, prefix='train')
                self._print_progress(epoch, train_info)
                # Avoid double-appending when the last epoch coincides with log_interval.
                if not history['epoch'] or history['epoch'][-1] != epoch:
                    history['epoch'].append(epoch)
                    for key in history_keys:
                        history[key].append(train_info.get(key, 0.0))
            
            # Validation
            run_validation = self.config['training'].get('run_validation', False)
            if run_validation and epoch % self.config['training']['validation_interval'] == 0:
                val_info = self.validate()
                self._log_metrics(val_info, prefix='val')
                lr = self.optimizer.param_groups[0]['lr']
                if self.manifold_type == 'sphere_product':
                    nij_str = '  '.join(
                        f"Nij{ia}{ib}: {val_info.get(f'nij_{ia}{ib}', 0.0):.4e}"
                        for ia, ib in [(0,0),(0,1),(1,0),(1,1)]
                    )
                    rank_str = ''
                    if self.config['logging'].get('compute_nijenhuis_rank', False):
                        rank_vals = [val_info.get(f'nij_rank_{ia}{ib}', float('nan'))
                                     for ia, ib in [(0,0),(0,1),(1,0),(1,1)]]
                        rank_str = ' | Rank: ' + '/'.join(f'{r:.1f}' for r in rank_vals)
                    print(f"  Val     {epoch:5d} | Loss: {val_info['total_loss']:.4e} | "
                          f"{nij_str} | Overlap: {val_info['overlap_norm']:.4e}{rank_str}")
                else:
                    rank_str = ''
                    if self.config['logging'].get('compute_nijenhuis_rank', False):
                        r1 = val_info.get('nij_rank1', float('nan'))
                        r2 = val_info.get('nij_rank2', float('nan'))
                        rank_str = f' | Rank: {r1:.1f}/{r2:.1f}'
                    print(f"  Val     {epoch:5d} | Loss: {val_info['total_loss']:.4e} | "
                          f"Nij1: {val_info['nij1']:.4e}  Nij2: {val_info['nij2']:.4e} | "
                          f"Overlap: {val_info['overlap_norm']:.4e}{rank_str}")

                # Check for improvement
                if val_info['total_loss'] < self.best_loss:
                    self.best_loss = val_info['total_loss']
                    self.patience_counter = 0
                    self._save_checkpoint('best_model.pt')
                else:
                    self.patience_counter += 1

                # Early stopping
                if (self.config['advanced']['early_stopping'] and
                        self.patience_counter >= self.config['advanced']['early_stopping_patience']):
                    print(f"\nEarly stopping at epoch {epoch}")
                    break
            elif not run_validation:
                # Track best train loss when validation is disabled.
                if train_info['total_loss'] < self.best_loss:
                    self.best_loss = train_info['total_loss']
                    self.patience_counter = 0
                    self._save_checkpoint('best_model.pt')
                else:
                    self.patience_counter += 1

                if (self.config['advanced']['early_stopping'] and
                        self.patience_counter >= self.config['advanced']['early_stopping_patience']):
                    print(f"\nEarly stopping at epoch {epoch}")
                    break
            
            # Checkpointing
            if epoch % checkpoint_interval == 0 and epoch > 0:
                self._save_checkpoint(f'checkpoint_epoch_{epoch}.pt')
            
            # Update scheduler
            if self.scheduler is not None and not using_lbfgs:
                self.scheduler.step()
        
        # Final checkpoint
        self._save_checkpoint('final_model.pt')
        print(f"\nTraining complete! Best loss: {self.best_loss:.6e}")

        self._evaluate_test()

        loss_json_path = self.output_dir / 'loss_data.json'
        with open(loss_json_path, 'w') as f:
            json.dump(self.loss_data, f, indent=2)
        print(f"Loss data saved to: {loss_json_path}")

        self._plot_loss_curves(history)
    
    @staticmethod
    def _block_diagonality(J_dict: dict, dim_alpha: int) -> float:
        """Compute off-diagonal-ness measure for sphere product J.

        For J = [[A, B], [C, D]] (n = dim_alpha, m = total_dim - n):

            ρ = (mean|B| + mean|C|) / (mean|A| + mean|D|)

        Averaged across all patches in J_dict.  Returns 0.0 when the
        denominator is zero (e.g. J = 0 everywhere).
        """
        n = dim_alpha
        ratios: list[float] = []
        for J in J_dict.values():
            A = J[:, :n, :n]
            B = J[:, :n, n:]
            C = J[:, n:, :n]
            D = J[:, n:, n:]
            num = B.abs().mean() + C.abs().mean()
            den = A.abs().mean() + D.abs().mean()
            ratios.append((num / den).item() if den.item() > 0.0 else 0.0)
        return float(sum(ratios) / len(ratios)) if ratios else 0.0

    def _evaluate_test(self):
        """Evaluate the model on a held-out test set and store results in loss_data."""
        test_size = self.config['logging'].get('test_size', 10000)
        self.model.eval()

        with torch.no_grad():
            if self.manifold_type == 'sphere_product':
                xa, xb = self.sphere_product.sample_points(test_size)
                xa, xb = xa.to(device=self.device, dtype=self.dtype), xb.to(device=self.device, dtype=self.dtype)
                J_dict, xa, xb, xa2, xb2 = self.model(xa, xb)
                _, info = self.loss_fn(J_dict, xa, xb, xa2, xb2, self.model,
                                       hard_cutoff=True, report_all=True)
                if self.config['logging'].get('compute_block_diagonality', False):
                    info['block_diagonality'] = self._block_diagonality(
                        J_dict, self.config['geometry']['real_dimension_alpha']
                    )
            else:
                x = self.sphere.sample_points_patch1(test_size).to(device=self.device, dtype=self.dtype)
                J1, J2, x1, x2 = self.model(x)
                _, info = self.loss_fn(J1, J2, x1, x2, self.model, hard_cutoff=True, report_all=True)

        if self.manifold_type == 'sphere_product':
            c1 = compute_first_chern_number(
                self.model,
                manifold_type='sphere_product',
                dim_alpha=self.config['geometry']['real_dimension_alpha'],
                dim_beta=self.config['geometry']['real_dimension_beta'],
                x_alpha=xa,
                x_beta=xb,
            )
        else:
            c1 = compute_first_chern_number(
                self.model,
                manifold_type='sphere',
                manifold_dim=self.sphere.dim,
                x=x,
            )

        self.model.train()

        result = {k: float(v) for k, v in info.items()}
        result['c1'] = float(c1)
        result['test_size'] = test_size
        self.loss_data['test'] = result

        print(f"\nTest evaluation ({test_size} samples):")
        print(f"  Total loss : {result['total_loss']:.6e}")
        if self.manifold_type == 'sphere_product':
            for key in ['nij_00', 'nij_01', 'nij_10', 'nij_11']:
                print(f"  {key}      : {result.get(key, 0.0):.6e}")
            if self.config['logging'].get('compute_nijenhuis_rank', False):
                for ia, ib in [(0,0),(0,1),(1,0),(1,1)]:
                    key = f'nij_rank_{ia}{ib}'
                    if key in result:
                        print(f"  rank {ia}{ib}    : {result[key]:.2f}")
        else:
            print(f"  Nij1       : {result['nij1']:.6e}")
            print(f"  Nij2       : {result['nij2']:.6e}")
            if self.config['logging'].get('compute_nijenhuis_rank', False):
                if 'nij_rank1' in result:
                    print(f"  Rank (p1)  : {result['nij_rank1']:.2f}")
                    print(f"  Rank (p2)  : {result['nij_rank2']:.2f}")
        print(f"  Overlap    : {result['overlap_norm']:.6e}")
        if 'block_diagonality' in result:
            print(f"  Block-diag ρ: {result['block_diagonality']:.6e}")
        print(f"  c₁          : {result['c1']:.6f}")

    def _log_metrics(self, info: dict, prefix: str):
        """Append metrics to loss_data."""
        record: dict = {'epoch': self.epoch}
        record.update({k: float(v) for k, v in info.items()})
        if prefix == 'train':
            record['lr'] = self.optimizer.param_groups[0]['lr']
        self.loss_data[prefix].append(record)
    
    def _print_progress(self, epoch: int, info: dict):
        """Print training progress."""
        lr = self.optimizer.param_groups[0]['lr']
        if self.manifold_type == 'sphere_product':
            nij_str = '  '.join(
                f"Nij{ia}{ib}: {info.get(f'nij_{ia}{ib}', 0.0):.4e}"
                for ia, ib in [(0,0),(0,1),(1,0),(1,1)]
            )
            print(f"Epoch {epoch:5d} | Loss: {info['total_loss']:.4e} | "
                  f"{nij_str} | Overlap: {info['overlap_norm']:.4e} | LR: {lr:.3e}")
        else:
            print(f"Epoch {epoch:5d} | Loss: {info['total_loss']:.4e} | "
                  f"Nij1: {info['nij1']:.4e}  Nij2: {info['nij2']:.4e} | "
                  f"Overlap: {info['overlap_norm']:.4e} | LR: {lr:.3e}")
    
    def _save_checkpoint(self, filename: str):
        """Save model checkpoint."""
        checkpoint = {
            'epoch': self.epoch,
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'optimizer_name': self.optimizer.__class__.__name__.lower(),
            'best_loss': self.best_loss,
            'config': self.config
        }
        
        if self.scheduler is not None:
            checkpoint['scheduler_state_dict'] = self.scheduler.state_dict()
        
        path = self.checkpoint_dir / filename
        torch.save(checkpoint, path)
    
    def _plot_loss_curves(self, history: dict):
        """Save loss curves to loss_curves.png after training."""
        epochs = history['epoch']
        weights = self.config.get('loss_weights', {})
        w_nij     = float(weights.get('nijenhuis', 1.0))
        w_overlap = float(weights.get('overlap', 1.0))

        plot_epochs = [e + 1 for e in epochs]

        fig, ax = plt.subplots(figsize=(9, 5))

        if self.manifold_type == 'sphere_product':
            nij_keys = ['nij_00', 'nij_01', 'nij_10', 'nij_11']
            nij_labels = ['Nij (00)', 'Nij (01)', 'Nij (10)', 'Nij (11)']
        else:
            nij_keys = ['nij1', 'nij2']
            nij_labels = ['Nij (patch 1)', 'Nij (patch 2)']

        active_nij   = w_nij     != 0.0 and all(k in history for k in nij_keys)
        active_ov    = w_overlap  != 0.0 and 'overlap_norm' in history

        if active_nij:
            for k, label in zip(nij_keys, nij_labels):
                ax.semilogy(plot_epochs, history[k], label=label, linestyle='--')
        if active_ov:
            ax.semilogy(plot_epochs, history['overlap_norm'], label='Overlap', linestyle=':')

        ax.set_xlabel('Epoch')
        ax.set_ylabel('Loss')
        ax.legend(loc='upper right')
        ax.grid(True, which='both', alpha=0.3)
        plt.tight_layout()
        out = self.output_dir / 'loss_curves.png'
        plt.savefig(out, dpi=150)
        plt.close()
        print(f"Loss curves saved to: {out}")

        if (
            self.manifold_type == 'sphere_product'
            and self.config['logging'].get('compute_block_diagonality', False)
            and 'block_diagonality' in history
        ):
            fig2, ax2 = plt.subplots(figsize=(8, 4))
            ax2.semilogy(plot_epochs, history['block_diagonality'],
                         label='Off-diagonal ratio \u03c1', color='C4')
            ax2.set_xlabel('Epoch')
            ax2.set_ylabel('\u03c1')
            ax2.legend(loc='upper right')
            ax2.grid(True, which='both', alpha=0.3)
            plt.tight_layout()
            out2 = self.output_dir / 'block_diagonality.png'
            plt.savefig(out2, dpi=150)
            plt.close()
            print(f"Block diagonality plot saved to: {out2}")

        if self.config['logging'].get('compute_nijenhuis_rank', False):
            val_records = self.loss_data.get('val', [])
            if self.manifold_type == 'sphere_product':
                rank_keys   = [f'nij_rank_{ia}{ib}' for ia, ib in [(0,0),(0,1),(1,0),(1,1)]]
                rank_labels = [f'Rank N ({ia}{ib})' for ia, ib in [(0,0),(0,1),(1,0),(1,1)]]
            else:
                rank_keys   = ['nij_rank1', 'nij_rank2']
                rank_labels = ['Rank N (patch 1)', 'Rank N (patch 2)']
            val_epochs = [r['epoch'] + 1 for r in val_records if rank_keys[0] in r]
            if val_epochs:
                fig3, ax3 = plt.subplots(figsize=(8, 4))
                for key, label in zip(rank_keys, rank_labels):
                    vals = [r[key] for r in val_records if key in r]
                    ax3.plot(val_epochs[:len(vals)], vals, label=label)
                ax3.set_xlabel('Epoch')
                ax3.set_ylabel('Mean numerical rank of N_J')
                ax3.legend(loc='upper right')
                ax3.grid(True, alpha=0.3)
                plt.tight_layout()
                out3 = self.output_dir / 'nijenhuis_rank.png'
                plt.savefig(out3, dpi=150)
                plt.close()
                print(f"Nijenhuis rank plot saved to: {out3}")

    def load_checkpoint(self, path: str):
        """Load model from checkpoint."""
        checkpoint = torch.load(path, map_location=self.device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        if checkpoint.get('optimizer_name') == 'lbfgs' and not isinstance(self.optimizer, optim.LBFGS):
            self.optimizer = self._create_lbfgs_optimizer()
            self.scheduler = None
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        self.epoch = checkpoint['epoch']
        self.best_loss = checkpoint['best_loss']
        
        if 'scheduler_state_dict' in checkpoint and self.scheduler is not None:
            self.scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
        
        print(f"Loaded checkpoint from epoch {self.epoch}")


def main():
    parser = argparse.ArgumentParser(
        description='Train PINN for almost complex structure'
    )
    parser.add_argument(
        '--config',
        type=str,
        default='configs/config_sphere2.yaml',
        help='Path to configuration file'
    )
    parser.add_argument(
        '--checkpoint',
        type=str,
        default=None,
        help='Path to checkpoint to resume from'
    )
    
    args = parser.parse_args()
    
    # Create trainer
    trainer = Trainer(args.config)
    
    # Load checkpoint if provided
    if args.checkpoint is not None:
        trainer.load_checkpoint(args.checkpoint)
    
    # Train
    trainer.train()


if __name__ == '__main__':
    main()
