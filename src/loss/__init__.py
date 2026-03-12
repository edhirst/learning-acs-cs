"""Loss package for Nijenhuis tensor and training losses."""

from .sphere_loss import (
    NijenhuisTensorLoss,
    OverlapLoss,
    TotalLoss,
)
from .product_loss import (
    ProductNijenhuisLoss,
    ProductOverlapLoss,
    ProductTotalLoss,
)

__all__ = [
    'NijenhuisTensorLoss',
    'OverlapLoss',
    'TotalLoss',
    'ProductNijenhuisLoss',
    'ProductOverlapLoss',
    'ProductTotalLoss',
]
