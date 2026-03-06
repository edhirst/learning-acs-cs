"""
Main training script for learning almost complex structures.
"""

import torch
import torch.optim as optim
import yaml
import json
import argparse
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

from src.geometry import EvenSphere
from src.models import TwoPatchAlmostComplexNet
from src.loss import TotalLoss


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
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(f"Using device: {self.device}")
        
        # Setup geometry
        self.sphere = EvenSphere(
            dimension=self.config['geometry']['dimension'],
            device=self.device
        )
        print(f"Manifold: S^{self.sphere.dim} (complex dimension {self.sphere.complex_dim})")
        
        # Setup model
        self.model = TwoPatchAlmostComplexNet(
            manifold_dim=self.sphere.dim,
            hidden_layers=self.config['network']['hidden_layers'],
            activation=self.config['network']['activation'],
            use_fourier_features=self.config['network']['use_fourier_features'],
            num_fourier_features=self.config['network']['num_fourier_features'],
            fourier_scale=self.config['network']['fourier_scale'],
            initialization=self.config['network']['initialization'],
        ).to(self.device)
        
        print(f"Model parameters: {sum(p.numel() for p in self.model.parameters())}")
        
        # Setup loss function
        self.loss_fn = TotalLoss(
            manifold_dim=self.sphere.dim,
            weights=self.config['loss_weights']
        ).to(self.device)
        
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
            return optim.lr_scheduler.CosineAnnealingWarmRestarts(
                self.optimizer,
                T_0=params['T_0'],
                T_mult=params.get('T_mult', 2),
                eta_min=params['eta_min']
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
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_dir = Path(self.config['logging']['save_dir']) / timestamp
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.checkpoint_dir = self.output_dir / 'checkpoints'
        self.checkpoint_dir.mkdir(exist_ok=True)

        self.loss_data: dict[str, list | dict] = {'train': [], 'val': [], 'test': {}}

        # Save config
        with open(self.output_dir / 'config.yaml', 'w') as f:
            yaml.dump(self.config, f)

        print(f"Outputs will be saved to: {self.output_dir}")
    
    def sample_batch(self, batch_size: int) -> torch.Tensor:
        """Sample a batch of patch-1 ball coordinates."""
        return self.sphere.sample_points_patch1(batch_size).to(self.device)
    
    def train_step(self, x: torch.Tensor) -> dict:
        """Perform one training step."""
        self.model.train()
        self.optimizer.zero_grad()

        # Forward pass
        J1, J2, x1, x2 = self.model(x)

        # Compute loss
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

        # Nijenhuis loss requires autograd through x; do not use no_grad here.
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
        dim = self.sphere.dim

        print(f"\nPretraining: A→0 warm-start ({pretrain_epochs} epochs, "
              f"{pretrain_points} points/step)")

        for subnet_name, subnet in [('patch1', self.model.patch1),
                                     ('patch2', self.model.patch2)]:
            opt = torch.optim.Adam(
                subnet.parameters(),
                lr=self.config['training']['learning_rate']
            )
            for epoch in range(pretrain_epochs):
                subnet.train()
                opt.zero_grad()
                x = self.sample_batch(pretrain_points)
                A = subnet.forward_raw(x)              # (N, 2n, 2n)
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
        history: dict[str, list] = {'epoch': [], 'total_loss': [], 'nij1': [], 'nij2': [], 'overlap_norm': []}
        
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
                    for key in ('total_loss', 'nij1', 'nij2', 'overlap_norm'):
                        history[key].append(train_info[key])
            
            # Validation
            run_validation = self.config['training'].get('run_validation', False)
            if run_validation and epoch % self.config['training']['validation_interval'] == 0:
                val_info = self.validate()
                self._log_metrics(val_info, prefix='val')
                lr = self.optimizer.param_groups[0]['lr']
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
        x = self.sphere.sample_points_patch1(test_size).to(self.device)
        J1, J2, x1, x2 = self.model(x)
        _, info = self.loss_fn(J1, J2, x1, x2, self.model, hard_cutoff=True, report_all=True)
        self.model.train()

        result = {k: float(v) for k, v in info.items()}
        result['test_size'] = test_size
        self.loss_data['test'] = result

        print(f"\nTest evaluation ({test_size} samples):")
        print(f"  Total loss : {result['total_loss']:.6e}")
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

        active = []
        if w_nij != 0.0:
            active += ['nij1', 'nij2']
        if w_overlap != 0.0:
            active += ['overlap_norm']

        # Shift epoch index to 1-based so x-axis runs 1 … N_epochs.
        plot_epochs = [e + 1 for e in epochs]

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.semilogy(plot_epochs, history['total_loss'], label='Total', linewidth=2)

        # Only add component lines when there is more than one active component
        # to avoid redundant overlapping lines.
        if len(active) > 1:
            if w_nij != 0.0:
                ax.semilogy(plot_epochs, history['nij1'],         label='Nij (patch 1)', linestyle='--')
                ax.semilogy(plot_epochs, history['nij2'],         label='Nij (patch 2)', linestyle='--')
            if w_overlap != 0.0:
                ax.semilogy(plot_epochs, history['overlap_norm'], label='Overlap',        linestyle=':')

        ax.set_xlabel('Epoch')
        ax.set_ylabel('Loss')
        ax.legend()
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
