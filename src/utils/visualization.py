"""
Utility functions for visualization and analysis.
"""

import torch
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
from typing import Optional, Tuple


def plot_constraint_violations(J: torch.Tensor, 
                               save_path: Optional[str] = None,
                               show: bool = True) -> plt.Figure:
    """
    Plot histograms of constraint violations.
    
    Args:
        J: Almost complex structure of shape (batch, dim, dim)
        save_path: Optional path to save figure
        show: Whether to display the figure
        
    Returns:
        matplotlib Figure
    """
    batch_size, dim, _ = J.shape
    
    # Compute constraint violations
    J_squared = torch.bmm(J, J)
    I = torch.eye(dim, device=J.device).unsqueeze(0).expand(batch_size, -1, -1)
    
    structure_error = torch.norm((J_squared + I).reshape(batch_size, -1), dim=1)
    
    # Eigenvalue magnitudes (should be 1 for ±i)
    eigenvalues = []
    for i in range(min(batch_size, 200)):
        eigvals = torch.linalg.eigvals(J[i].cpu())
        eigenvalues.append(torch.abs(eigvals))
    eigmags = torch.stack(eigenvalues).flatten().numpy()
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Structure constraint
    axes[0].hist(structure_error.cpu().numpy(), bins=50, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('||J² + I||', fontsize=12)
    axes[0].set_ylabel('Frequency', fontsize=12)
    axes[0].set_title('Structure Constraint Violation', fontsize=14)
    axes[0].grid(alpha=0.3)
    
    # Add statistics
    mean_str = structure_error.mean().item()
    max_str = structure_error.max().item()
    axes[0].axvline(mean_str, color='r', linestyle='--', label=f'Mean: {mean_str:.2e}')
    axes[0].legend()
    
    # Eigenvalue magnitude distribution
    axes[1].hist(eigmags, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    axes[1].axvline(1.0, color='r', linestyle='--', label='Expected (|±i| = 1)')
    axes[1].set_xlabel('|λ|', fontsize=12)
    axes[1].set_ylabel('Frequency', fontsize=12)
    axes[1].set_title('Eigenvalue Magnitude Distribution', fontsize=14)
    axes[1].grid(alpha=0.3)
    axes[1].legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    if show:
        plt.show()
    
    return fig


def plot_sphere_with_errors(x: torch.Tensor,
                            J: torch.Tensor,
                            save_path: Optional[str] = None,
                            show: bool = True) -> plt.Figure:
    """
    Plot sphere with constraint errors (only for S²).
    
    Args:
        x: Points on S² of shape (batch, 3)
        J: Almost complex structure of shape (batch, 2, 2)
        save_path: Optional path to save figure
        show: Whether to display the figure
        
    Returns:
        matplotlib Figure
    """
    if x.shape[1] != 3:
        raise ValueError("This function only works for S² (3D embedding)")
    
    batch_size = x.shape[0]
    
    # Compute constraint errors
    J_squared = torch.bmm(J, J)
    I = torch.eye(2, device=J.device).unsqueeze(0).expand(batch_size, -1, -1)
    errors = torch.norm((J_squared + I).reshape(batch_size, -1), dim=1)
    
    # Create 3D plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    x_np = x.cpu().numpy()
    errors_np = errors.cpu().numpy()
    
    scatter = ax.scatter(x_np[:, 0], x_np[:, 1], x_np[:, 2],
                        c=errors_np, cmap='viridis', s=20)
    
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel('Z', fontsize=12)
    ax.set_title('S² with Constraint Errors', fontsize=14)
    
    cbar = plt.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label('||J² + I||', fontsize=12)
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    if show:
        plt.show()
    
    return fig


def plot_eigenvalue_distribution(J: torch.Tensor,
                                 save_path: Optional[str] = None,
                                 show: bool = True) -> plt.Figure:
    """
    Plot eigenvalue distribution in complex plane.
    
    For a perfect almost complex structure, eigenvalues should be ±i.
    
    Args:
        J: Almost complex structure of shape (batch, dim, dim)
        save_path: Optional path to save figure
        show: Whether to display the figure
        
    Returns:
        matplotlib Figure
    """
    batch_size = J.shape[0]
    
    # Compute eigenvalues
    eigenvalues = []
    for i in range(min(batch_size, 500)):  # Limit for efficiency
        eigvals = torch.linalg.eigvals(J[i].cpu())
        eigenvalues.append(eigvals)
    
    eigenvalues = torch.stack(eigenvalues).numpy()
    
    # Create plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Complex plane scatter
    axes[0].scatter(eigenvalues.real, eigenvalues.imag, alpha=0.3, s=10)
    axes[0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axes[0].axvline(x=0, color='k', linestyle='--', alpha=0.3)
    
    # Unit circle
    theta = np.linspace(0, 2*np.pi, 100)
    axes[0].plot(np.cos(theta), np.sin(theta), 'r--', alpha=0.5, label='Unit circle')
    
    # Expected eigenvalues ±i
    axes[0].plot([0], [1], 'ro', markersize=10, label='Expected (+i)')
    axes[0].plot([0], [-1], 'ro', markersize=10, label='Expected (-i)')
    
    axes[0].set_xlabel('Real', fontsize=12)
    axes[0].set_ylabel('Imaginary', fontsize=12)
    axes[0].set_title('Eigenvalues in Complex Plane', fontsize=14)
    axes[0].legend()
    axes[0].axis('equal')
    axes[0].grid(alpha=0.3)
    
    # Magnitude histogram
    magnitudes = np.abs(eigenvalues).flatten()
    axes[1].hist(magnitudes, bins=50, edgecolor='black', alpha=0.7)
    axes[1].axvline(1.0, color='r', linestyle='--', label='Expected (|±i| = 1)')
    axes[1].set_xlabel('|λ|', fontsize=12)
    axes[1].set_ylabel('Frequency', fontsize=12)
    axes[1].set_title('Eigenvalue Magnitude Distribution', fontsize=14)
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    if show:
        plt.show()
    
    return fig


def plot_training_history(log_file: str,
                          save_path: Optional[str] = None,
                          show: bool = True) -> plt.Figure:
    """
    Plot training history from log file.
    
    Args:
        log_file: Path to training log
        save_path: Optional path to save figure
        show: Whether to display the figure
        
    Returns:
        matplotlib Figure
    """
    # This is a placeholder - in practice, you'd parse TensorBoard logs
    # or save training metrics to a CSV
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Example plot structure
    axes[0, 0].set_title('Total Loss')
    axes[0, 1].set_title('Nijenhuis Norm')
    axes[1, 0].set_title('Structure Constraint')
    axes[1, 1].set_title('Smoothness')
    
    for ax in axes.flat:
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Loss')
        ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    if show:
        plt.show()
    
    return fig


def compute_metrics_summary(J: torch.Tensor) -> dict:
    """
    Compute summary statistics for almost complex structure.
    
    Args:
        J: Almost complex structure of shape (batch, dim, dim)
        
    Returns:
        Dictionary of metrics
    """
    batch_size, dim, _ = J.shape
    
    # Structure constraint
    J_squared = torch.bmm(J, J)
    I = torch.eye(dim, device=J.device).unsqueeze(0).expand(batch_size, -1, -1)
    structure_error = torch.norm((J_squared + I).reshape(batch_size, -1), dim=1)
    

    # Frobenius norms
    j_norm = torch.norm(J.reshape(batch_size, -1), dim=1)
    
    metrics = {
        'structure_error_mean': structure_error.mean().item(),
        'structure_error_std': structure_error.std().item(),
        'structure_error_max': structure_error.max().item(),
        'structure_error_min': structure_error.min().item(),
        'j_norm_mean': j_norm.mean().item(),
        'j_norm_std': j_norm.std().item(),
    }
    
    return metrics


def print_metrics_table(metrics: dict):
    """Print metrics in a formatted table."""
    print("\n" + "="*60)
    print("METRICS SUMMARY")
    print("="*60)
    
    print("\nStructure Constraint ||J² + I||:")
    print(f"  Mean: {metrics['structure_error_mean']:.6e}")
    print(f"  Std:  {metrics['structure_error_std']:.6e}")
    print(f"  Max:  {metrics['structure_error_max']:.6e}")
    print(f"  Min:  {metrics['structure_error_min']:.6e}")
    
    print("\nJ Frobenius Norm:")
    print(f"  Mean: {metrics['j_norm_mean']:.6e}")
    print(f"  Std:  {metrics['j_norm_std']:.6e}")
    
    print("="*60 + "\n")


def plot_chern_curvature(
    x: torch.Tensor,
    chern_density: torch.Tensor,
    nij_norm: Optional[torch.Tensor] = None,
    save_path: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """
    Plot Chern curvature density in ball coordinates.

    For S² (2D ball coords), produces a 2D scatter plot of the Chern density
    (1/2π) F_{01}(x) √det(g(x)) across the unit disc.  Optionally overlays
    the Nijenhuis norm magnitude as point sizes.

    For higher-dimensional manifolds, plots the marginal distribution over the
    first two ball coordinates, coloured by Chern density.

    Args:
        x:             Ball coordinates, shape (N, d).
        chern_density: Chern curvature density at each point, shape (N,).
        nij_norm:      Optional Nijenhuis norm at each point, shape (N,).
                       When given, used as scatter point sizes.
        save_path:     Optional path to save figure.
        show:          Whether to call plt.show().

    Returns:
        matplotlib Figure.
    """
    x_np = x.detach().cpu().numpy()
    c_np = chern_density.detach().cpu().numpy()

    ncols = 2 if nij_norm is not None else 1
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5))
    if ncols == 1:
        axes = [axes]

    # --- Chern density scatter ---
    ax = axes[0]
    vmax = np.abs(c_np).max() + 1e-12
    sc = ax.scatter(
        x_np[:, 0], x_np[:, 1],
        c=c_np, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
        s=6, alpha=0.7,
    )
    plt.colorbar(sc, ax=ax, label='Chern density')
    # Draw unit circle boundary
    theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(theta), np.sin(theta), 'k--', lw=0.8, alpha=0.4)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect('equal')
    ax.set_xlabel('x₀')
    ax.set_ylabel('x₁')
    ax.set_title('Chern curvature density (1/2π) F₀₁ √det(g)')
    ax.grid(alpha=0.2)

    if nij_norm is not None:
        # --- Nijenhuis norm scatter ---
        nij_np = nij_norm.detach().cpu().numpy()
        ax2 = axes[1]
        sc2 = ax2.scatter(
            x_np[:, 0], x_np[:, 1],
            c=nij_np, cmap='viridis',
            s=6, alpha=0.7,
        )
        plt.colorbar(sc2, ax=ax2, label='Nijenhuis norm')
        ax2.plot(np.cos(theta), np.sin(theta), 'k--', lw=0.8, alpha=0.4)
        ax2.set_xlim(-1.05, 1.05)
        ax2.set_ylim(-1.05, 1.05)
        ax2.set_aspect('equal')
        ax2.set_xlabel('x₀')
        ax2.set_ylabel('x₁')
        ax2.set_title('Nijenhuis norm ‖N_J‖')
        ax2.grid(alpha=0.2)

    fig.tight_layout()

    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig


if __name__ == "__main__":
    # Example usage
    print("Testing visualization utilities...")
    
    # Create random almost complex structure
    batch_size = 100
    dim = 4
    
    J = torch.randn(batch_size, dim, dim)
    J = (J - J.transpose(1, 2)) / 2  # Make antisymmetric
    
    # Compute metrics
    metrics = compute_metrics_summary(J)
    print_metrics_table(metrics)
    
    # Create visualizations
    print("Creating visualizations...")
    plot_constraint_violations(J, show=False)
    plot_eigenvalue_distribution(J, show=False)
    
    print("✓ Visualization utilities working correctly!")
