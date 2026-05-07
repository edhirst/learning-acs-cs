"""Geometry package for manifold operations."""

from .sphere import EvenSphere, RADIAL_MIDPOINT, sphere_uniform_sample
from .sphere_product import SphereProduct
from .hopf_complex_structure import hopf_complex_structure, hopf_j_all_patches

__all__ = [
    'EvenSphere',
    'RADIAL_MIDPOINT',
    'sphere_uniform_sample',
    'SphereProduct',
    'hopf_complex_structure',
    'hopf_j_all_patches',
]
