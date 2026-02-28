"""Utils package for visualization and utilities."""

from .visualization import (
    plot_constraint_violations,
    plot_sphere_with_errors,
    plot_eigenvalue_distribution,
    compute_metrics_summary,
    print_metrics_table
)

__all__ = [
    'plot_constraint_violations',
    'plot_sphere_with_errors',
    'plot_eigenvalue_distribution',
    'compute_metrics_summary',
    'print_metrics_table'
]
