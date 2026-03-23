"""
Main training script for learning almost complex structures.
"""

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
from src.models import TwoPatchAlmostComplexNet, FourPatchProductNet
from src.loss import TotalLoss
from src.loss.product_loss import ProductTotalLoss


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

        # Setup dtype
        dtype_str = self.config['training'].get('dtype', 'float64')
        self.dtype = torch.float64 if dtype_str == 'float64' else torch.float32
        if self.dtype == torch.float64 and self.device.type == 'mps':
            print("Warning: MPS does not support float64; falling back to CPU.")
            self.device = torch.device('cpu')
        print(f"Using dtype : {dtype_str}")
        
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
            ).to(self.device).to(self.dtype)

            self.loss_fn = ProductTotalLoss(
                dim_alpha=dim_alpha,
                dim_beta=dim_beta,
                weights=self.config['loss_weights'],
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
            ).to(self.device).to(self.dtype)

            self.loss_fn = TotalLoss(
                manifold_dim=self.sphere.dim,
                weights=self.config['loss_weights'],
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
    
    def train_step(self, x) -> dict:
        """Perform one training step."""
        self.model.train()
        self.optimizer.zero_grad()

        if self.manifold_type == 'sphere_product':
            x_alpha, x_beta = x
            J_dict, x_alpha, x_beta, x_a2, x_b2 = self.model(x_alpha, x_beta)
            loss, info = self.loss_fn(J_dict, x_alpha, x_beta, x_a2, x_b2, self.model)
        else:
            J1, J2, x1, x2 = self.model(x)
            loss, info = self.loss_fn(J1, J2, x1, x2, self.model)
        
        # Backward pass
        loss.backward()
        
        # Gradient clipping
        if self.config['advanced']['gradient_clipping'] > 0:
            torch.nn.utils.clip_grad_norm_(
                self.model.parameters(),
                self.config['advanced']['gradient_clipping']
            )
        
        self.optimizer.step()
        
        return info
    
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

    def train(self):
        """Main training loop."""
        self.pretrain()

        print("\nStarting training...")
        print(f"Epochs: {self.config['training']['epochs']}")
        print(f"Batch size: {self.config['training']['batch_size']}")
        
        num_epochs = self.config['training']['epochs']
        batch_size = self.config['training']['batch_size']
        log_interval = self.config['logging']['log_interval']
        checkpoint_interval = self.config['logging']['checkpoint_interval']
        if self.manifold_type == 'sphere_product':
            history_keys = ['total_loss', 'nij_00', 'nij_01', 'nij_10', 'nij_11', 'overlap_norm']
        else:
            history_keys = ['total_loss', 'nij1', 'nij2', 'overlap_norm']
        history: dict[str, list] = {'epoch': [], **{k: [] for k in history_keys}}
        
        for epoch in range(num_epochs):
            self.epoch = epoch
            
            # Sample batch
            x = self.sample_batch(batch_size)
            
            # Training step
            train_info = self.train_step(x)
            
            # Logging
            is_last_epoch = (epoch == num_epochs - 1)
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
                    print(f"  Val     {epoch:5d} | Loss: {val_info['total_loss']:.4e} | "
                          f"{nij_str} | Overlap: {val_info['overlap_norm']:.4e}")
                else:
                    print(f"  Val     {epoch:5d} | Loss: {val_info['total_loss']:.4e} | "
                          f"Nij1: {val_info['nij1']:.4e}  Nij2: {val_info['nij2']:.4e} | "
                          f"Overlap: {val_info['overlap_norm']:.4e}")

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
            if self.scheduler is not None:
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
            else:
                x = self.sphere.sample_points_patch1(test_size).to(device=self.device, dtype=self.dtype)
                J1, J2, x1, x2 = self.model(x)
                _, info = self.loss_fn(J1, J2, x1, x2, self.model, hard_cutoff=True, report_all=True)

        self.model.train()

        result = {k: float(v) for k, v in info.items()}
        result['test_size'] = test_size
        self.loss_data['test'] = result

        print(f"\nTest evaluation ({test_size} samples):")
        print(f"  Total loss : {result['total_loss']:.6e}")
        if self.manifold_type == 'sphere_product':
            for key in ['nij_00', 'nij_01', 'nij_10', 'nij_11']:
                print(f"  {key}      : {result.get(key, 0.0):.6e}")
        else:
            print(f"  Nij1       : {result['nij1']:.6e}")
            print(f"  Nij2       : {result['nij2']:.6e}")
        print(f"  Overlap    : {result['overlap_norm']:.6e}")

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
    
    def load_checkpoint(self, path: str):
        """Load model from checkpoint."""
        checkpoint = torch.load(path, map_location=self.device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
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
