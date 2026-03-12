"""
Neural network architectures for learning almost complex structures.

Two manifold types are supported:

sphere
------
`TwoPatchAlmostComplexNet` wraps two `PatchSubNet` instances covering S^{2n}
via stereographic ball patches.  Given patch-1 coordinates x₁ ∈ B(0,1), it
computes x₂ = φ(x₁) and evaluates both sub-networks.

sphere_product
--------------
`FourPatchProductNet` covers S^{m_α} × S^{m_β} with four patches
(i_α, i_β) ∈ {0,1}², one `PatchSubNet` per patch combination.  Each subnet
receives the full concatenated ball coordinates x = (x_α, x_β).

Per-patch J₀ for the product:
    The sign flip per patch reflects the orientation of the equatorial Jacobian:
    det(F_eq) = (-1)^{dim_k + 1} for each factor k.  For even dim_k this equals
    -1 (orientation-reversing → J_k flips sign), for odd dim_k this equals +1
    (orientation-preserving → no sign flip).

    Both dims even → j0^{(ia,ib)} = block( (-1)^{ia} J₀^α, (-1)^{ib} J₀^β )
    Both dims odd  → j0^{(ia,ib)} = J₀(d)  for all (ia,ib)  [same for all patches]

In all cases the exp(A) J₀ exp(−A) parametrisation guarantees J² = −I exactly
for any network weights.
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Dict, List, Optional, Tuple

from src.geometry.sphere import EvenSphere


def _matrix_exp(A: torch.Tensor) -> torch.Tensor:
    """torch.matrix_exp with CPU fallback for MPS (not yet implemented there)."""
    if A.device.type == 'mps':
        return torch.matrix_exp(A.cpu()).to(A.device)
    return torch.matrix_exp(A)


# ─────────────────────────────────────────────────────────────────────────────
# Shared building blocks
# ─────────────────────────────────────────────────────────────────────────────

class SinActivation(nn.Module):
    """Sine activation function."""

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return torch.sin(x)


class FourierFeatures(nn.Module):
    """Random Fourier feature encoding for positional encoding."""

    def __init__(self, input_dim: int, num_features: int, scale: float = 1.0):
        super().__init__()
        B = torch.randn(input_dim, num_features) * scale
        self.register_buffer('B', B)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: Input of shape (batch, input_dim).

        Returns:
            Encoded features of shape (batch, 2 * num_features).
        """
        proj = 2 * np.pi * x @ self.B
        return torch.cat([torch.sin(proj), torch.cos(proj)], dim=-1)


class PatchSubNet(nn.Module):
    """
    Sub-network for a single ball patch.

    Takes ball-patch coordinates x ∈ B^{2n} ⊂ ℝ^{2n} as input and returns
    an almost complex structure J = exp(A(x)) J₀ exp(−A(x)) ∈ ℝ^{2n×2n}.

    J² = −I holds exactly for any network weights by construction.
    """

    def __init__(
        self,
        manifold_dim: int,
        hidden_layers: List[int],
        activation: str = 'tanh',
        use_fourier_features: bool = True,
        num_fourier_features: int = 64,
        fourier_scale: float = 1.0,
        initialization: str = 'xavier',
        j0_sign: int = 1,
        j0_matrix: Optional[torch.Tensor] = None,
    ):
        """
        Args:
            manifold_dim: Real dimension of the manifold (must be even); also the
                input coordinate dimension.
            hidden_layers: Sizes of hidden layers.
            activation: Activation function name ('tanh', 'relu', 'gelu', 'sin').
            use_fourier_features: Whether to use Fourier feature encoding.
            num_fourier_features: Number of Fourier basis functions.
            fourier_scale: Scale for random Fourier features.
            initialization: Weight initialisation ('xavier', 'kaiming', 'orthogonal').
            j0_sign: Sign of the canonical structure J₀ (+1 or −1).  Ignored when
                j0_matrix is provided.  Patch 1 uses +J₀ (north-pole chart
                orientation); patch 2 uses −1 because the north→south chart
                transition x → x/|x|² is complex-anti-linear, mapping J → −J.
            j0_matrix: Optional explicit J₀ ∈ ℝ^{d×d} satisfying J₀² = −I.
                When provided, overrides j0_sign and _build_J0.  Used by
                FourPatchProductNet to supply per-patch block-diagonal J₀s.
        """
        super().__init__()

        if manifold_dim % 2 != 0:
            raise ValueError("manifold_dim must be even")

        self.manifold_dim = manifold_dim
        self.complex_dim  = manifold_dim // 2
        n = manifold_dim // 2
        # Effective DOFs: anti-commuting subspace {A : A J₀ + J₀ A = 0} has
        # real dimension 2n² (half of gl(2n,R)).  The complementary commuting
        # component is a gauge freedom: exp(A) J₀ exp(−A) is unchanged when A
        # is shifted by any matrix that commutes with J₀.
        self.num_dof = 2 * n * n

        # Canonical complex structure: use the supplied matrix when given,
        # otherwise build the standard block-diagonal J₀ (with sign).
        if j0_matrix is not None:
            J0 = j0_matrix
        else:
            J0 = j0_sign * self._build_J0(manifold_dim)
        self.register_buffer('J0', J0)

        # Orthonormal basis for the anti-commuting subspace, shape (num_dof, dim²).
        # Each row is a flattened basis element; used to embed network outputs
        # into the gauge-fixed subspace of gl(2n,R).
        basis = self._build_anticommuting_basis(manifold_dim, J0)
        self.register_buffer('basis', basis)  # (num_dof, dim²)

        # Input encoding
        if use_fourier_features:
            self.fourier = FourierFeatures(manifold_dim, num_fourier_features, fourier_scale)
            input_dim = 2 * num_fourier_features
        else:
            self.fourier = None
            input_dim = manifold_dim

        # Activation
        act_map = {
            'tanh': nn.Tanh,
            'relu': nn.ReLU,
            'gelu': nn.GELU,
            'sin':  SinActivation,
        }
        if activation not in act_map:
            raise ValueError(f"Unknown activation: {activation}")
        act_cls = act_map[activation]

        # MLP layers
        layers = []
        prev = input_dim
        for h in hidden_layers:
            layers += [nn.Linear(prev, h), act_cls()]
            prev = h
        layers.append(nn.Linear(prev, self.num_dof))
        self.network = nn.Sequential(*layers)

        self._initialize_weights(initialization)

    @staticmethod
    def _build_anticommuting_basis(dim: int, J0: torch.Tensor) -> torch.Tensor:
        """
        Orthonormal basis for {A ∈ gl(2n,ℝ) : A J₀ + J₀ A = 0}.

        Builds the anti-commutator map M on gl(2n,ℝ) and extracts its null
        space via SVD.  The result has shape (2n², (2n)²); each row is a
        flattened basis element.

        Args:
            dim: Even real dimension 2n.
            J0:  Canonical complex structure, shape (dim, dim).

        Returns:
            basis of shape (2n², dim²), dtype float32.
        """
        J0_np = J0.numpy()
        d = dim
        # Build M column-by-column: column i*d+j is the flattened image of the
        # basis matrix e_{ij} under A ↦ A J₀ + J₀ A.
        M = np.zeros((d * d, d * d), dtype=np.float64)
        for i in range(d):
            for j in range(d):
                E = np.zeros((d, d), dtype=np.float64)
                E[i, j] = 1.0
                M[:, i * d + j] = (E @ J0_np + J0_np @ E).flatten()
        # Right null space of M = anti-commuting subspace.
        _, s, Vt = np.linalg.svd(M, full_matrices=True)
        null_basis = Vt[s < 1e-8]     # (2n², dim²)
        return torch.tensor(null_basis, dtype=torch.float32)

    @staticmethod
    def _build_J0(dim: int) -> torch.Tensor:
        """
        Build the canonical complex structure J₀ = block_diag([[0,−1],[1,0]], …).

        Satisfies J₀² = −I.

        Args:
            dim: Even real dimension.

        Returns:
            J₀ of shape (dim, dim).
        """
        J0 = torch.zeros(dim, dim)
        for k in range(dim // 2):
            i, j = 2 * k, 2 * k + 1
            J0[i, j] = -1.0
            J0[j, i] =  1.0
        return J0

    def _initialize_weights(self, method: str) -> None:
        for m in self.modules():
            if isinstance(m, nn.Linear):
                if method == 'xavier':
                    nn.init.xavier_normal_(m.weight)
                elif method == 'kaiming':
                    nn.init.kaiming_normal_(m.weight)
                elif method == 'orthogonal':
                    nn.init.orthogonal_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Compute almost complex structure J at ball-patch coordinates x.

        J = exp(A(x)) J₀ exp(−A(x))  satisfies J² = −I exactly.

        Args:
            x: Ball-patch coordinates, shape (batch, manifold_dim).

        Returns:
            J: Shape (batch, manifold_dim, manifold_dim).
        """
        features = self.fourier(x) if self.fourier is not None else x
        params = self.network(features)                      # (batch, num_dof)
        # Embed into the anti-commuting subspace of gl(2n,R).
        A_flat = params @ self.basis                         # (batch, dim²)
        batch = params.shape[0]
        A = A_flat.reshape(batch, self.manifold_dim, self.manifold_dim)

        expA  = _matrix_exp(A)
        expnA = _matrix_exp(-A)
        J0 = self.J0.unsqueeze(0).expand(batch, -1, -1)
        return expA @ J0 @ expnA

    def forward_raw(self, x: torch.Tensor) -> torch.Tensor:
        """Return A(x): the anti-commuting matrix before the Cayley-like transform.

        A lies exactly in {A : A J₀ + J₀ A = 0} by construction.

        Args:
            x: Ball-patch coordinates, shape (batch, manifold_dim).

        Returns:
            A: Shape (batch, manifold_dim, manifold_dim).
        """
        features = self.fourier(x) if self.fourier is not None else x
        params = self.network(features)                      # (batch, num_dof)
        A_flat = params @ self.basis                         # (batch, dim²)
        return A_flat.reshape(params.shape[0], self.manifold_dim, self.manifold_dim)


# ─────────────────────────────────────────────────────────────────────────────
# sphere manifold_type
# ─────────────────────────────────────────────────────────────────────────────

class TwoPatchAlmostComplexNet(nn.Module):
    """
    Global almost complex structure network with two stereographic ball patches.

    A single batch of patch-1 ball coordinates is given.  The network
    internally maps to patch-2 coordinates, evaluates both sub-networks,
    and returns the J matrices and coordinates for both patches.

    Together, the two patches cover all of S^{2n}.
    """

    def __init__(
        self,
        manifold_dim: int,
        hidden_layers: List[int],
        activation: str = 'tanh',
        use_fourier_features: bool = True,
        num_fourier_features: int = 64,
        fourier_scale: float = 1.0,
        initialization: str = 'xavier',
    ):
        """
        Args:
            manifold_dim: Real dimension of the sphere (2n).
            hidden_layers: Hidden layer sizes for each patch sub-network.
            activation: Activation function name.
            use_fourier_features: Whether to use Fourier features.
            num_fourier_features: Number of Fourier basis functions.
            fourier_scale: Scale for random Fourier features.
            initialization: Weight initialisation method.
        """
        super().__init__()

        subnet_kwargs = dict(
            manifold_dim=manifold_dim,
            hidden_layers=hidden_layers,
            activation=activation,
            use_fourier_features=use_fourier_features,
            num_fourier_features=num_fourier_features,
            fourier_scale=fourier_scale,
            initialization=initialization,
        )
        self.patch1 = PatchSubNet(**subnet_kwargs, j0_sign=+1)
        self.patch2 = PatchSubNet(**subnet_kwargs, j0_sign=-1)

        self.manifold_dim = manifold_dim
        self.complex_dim  = manifold_dim // 2
        self.num_dof      = manifold_dim * manifold_dim

    def forward(self, x1: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor,
                                                   torch.Tensor, torch.Tensor]:
        """
        Compute J in both patches given patch-1 ball coordinates.

        Args:
            x1: Patch-1 ball coordinates, shape (batch, manifold_dim).

        Returns:
            J1: J in patch 1, shape (batch, manifold_dim, manifold_dim).
            J2: J in patch 2, shape (batch, manifold_dim, manifold_dim).
            x1: Patch-1 coordinates (pass-through), shape (batch, manifold_dim).
            x2: Patch-2 coordinates, shape (batch, manifold_dim).
        """
        x2 = EvenSphere.patch_change_coords(x1)
        J1 = self.patch1(x1)
        J2 = self.patch2(x2)
        return J1, J2, x1, x2


# ─────────────────────────────────────────────────────────────────────────────
# sphere_product manifold_type
# ─────────────────────────────────────────────────────────────────────────────

# Canonical patch indices for product manifold
_PATCH_KEYS = [(0, 0), (0, 1), (1, 0), (1, 1)]


def _build_patch_j0s(
    dim_alpha: int,
    dim_beta: int,
) -> Dict[Tuple[int, int], torch.Tensor]:
    """
    Compute the per-patch canonical ACS J₀ for S^{m_α} × S^{m_β}.

    Args:
        dim_alpha: Real dimension of S^{m_α}.
        dim_beta:  Real dimension of S^{m_β}.

    Returns:
        Dict mapping (i_α, i_β) → J₀ tensor of shape (d, d), d = m_α + m_β.
    """
    d = dim_alpha + dim_beta
    j0s: Dict[Tuple[int, int], torch.Tensor] = {}

    if dim_alpha % 2 == 0 and dim_beta % 2 == 0:
        # Block-diagonal: j0^{(ia,ib)} = block( s_α J₀^α, s_β J₀^β )
        J0a = PatchSubNet._build_J0(dim_alpha)   # (m_α, m_α)
        J0b = PatchSubNet._build_J0(dim_beta)    # (m_β, m_β)
        for ia, ib in _PATCH_KEYS:
            J0 = torch.zeros(d, d)
            J0[:dim_alpha, :dim_alpha] = ((-1) ** ia) * J0a
            J0[dim_alpha:, dim_alpha:] = ((-1) ** ib) * J0b
            j0s[(ia, ib)] = J0
    else:
        # Both dims odd (even total).  No natural block decomposition; use
        # the standard J₀(d) on the full space.
        # For an odd-dimensional factor the patch-change Jacobian at the equator
        # has det = +1 (orientation-preserving), so F J F⁻¹ = J there — no sign
        # flip between patches.  All four patches share the same J₀.
        J0_full = PatchSubNet._build_J0(d)
        for ia, ib in _PATCH_KEYS:
            j0s[(ia, ib)] = J0_full

    return j0s


class FourPatchProductNet(nn.Module):
    """
    Global almost complex structure network for S^{m_α} × S^{m_β}.

    Contains four PatchSubNets, one per patch combination (i_α, i_β).
    Each subnet receives the full concatenated ball coordinates
    x = (x_α, x_β) ∈ ℝ^{m_α+m_β} as input.
    """

    def __init__(
        self,
        dim_alpha: int,
        dim_beta: int,
        hidden_layers: List[int],
        activation: str = 'tanh',
        use_fourier_features: bool = False,
        num_fourier_features: int = 64,
        fourier_scale: float = 1.0,
        initialization: str = 'xavier',
    ):
        """
        Args:
            dim_alpha: Real dimension of S^{m_α}.
            dim_beta:  Real dimension of S^{m_β}.
            hidden_layers: Hidden layer sizes shared by all four subnets.
            activation: Activation function ('tanh', 'relu', 'gelu', 'sin').
            use_fourier_features: Whether to use random Fourier features.
            num_fourier_features: Number of Fourier basis functions.
            fourier_scale: Scale for Fourier feature random matrix.
            initialization: Weight initialisation method.
        """
        super().__init__()

        self.dim_alpha  = dim_alpha
        self.dim_beta   = dim_beta
        self.total_dim  = dim_alpha + dim_beta
        self.complex_dim = self.total_dim // 2

        j0s = _build_patch_j0s(dim_alpha, dim_beta)

        subnet_kwargs = dict(
            manifold_dim=self.total_dim,
            hidden_layers=hidden_layers,
            activation=activation,
            use_fourier_features=use_fourier_features,
            num_fourier_features=num_fourier_features,
            fourier_scale=fourier_scale,
            initialization=initialization,
        )

        self.patches = nn.ModuleDict({
            f'patch_{ia}{ib}': PatchSubNet(**subnet_kwargs, j0_matrix=j0s[(ia, ib)])
            for ia, ib in _PATCH_KEYS
        })

    @property
    def subnets(self) -> List[PatchSubNet]:
        """All four PatchSubNets in canonical patch order."""
        return [self.patches[f'patch_{ia}{ib}'] for ia, ib in _PATCH_KEYS]

    def forward(
        self,
        x_alpha: torch.Tensor,
        x_beta:  torch.Tensor,
    ) -> Tuple[Dict[Tuple[int, int], torch.Tensor], torch.Tensor, torch.Tensor,
               torch.Tensor, torch.Tensor]:
        """
        Compute J in all four patches.

        Args:
            x_alpha: Patch-0 α-coords, shape (batch, dim_alpha).
            x_beta:  Patch-0 β-coords, shape (batch, dim_beta).

        Returns:
            J_dict: Dict (i_α, i_β) → J tensor of shape (batch, total_dim, total_dim).
            x_alpha: Pass-through, shape (batch, dim_alpha).
            x_beta:  Pass-through, shape (batch, dim_beta).
            x_a2:    φ_α(x_alpha), shape (batch, dim_alpha).
            x_b2:    φ_β(x_beta),  shape (batch, dim_beta).
        """
        x_a2 = EvenSphere.patch_change_coords(x_alpha)
        x_b2 = EvenSphere.patch_change_coords(x_beta)

        patch_alpha_coords = {0: x_alpha, 1: x_a2}
        patch_beta_coords  = {0: x_beta,  1: x_b2}

        J_dict: Dict[Tuple[int, int], torch.Tensor] = {}
        for ia, ib in _PATCH_KEYS:
            x_combined = torch.cat([patch_alpha_coords[ia], patch_beta_coords[ib]], dim=1)
            J_dict[(ia, ib)] = self.patches[f'patch_{ia}{ib}'](x_combined)

        return J_dict, x_alpha, x_beta, x_a2, x_b2

    def forward_raw(
        self,
        x_alpha: torch.Tensor,
        x_beta:  torch.Tensor,
    ) -> Dict[Tuple[int, int], torch.Tensor]:
        """
        Return A(x) for all four patches (used for pretraining, A → 0).

        Args:
            x_alpha: Patch-0 α-coords, shape (batch, dim_alpha).
            x_beta:  Patch-0 β-coords, shape (batch, dim_beta).

        Returns:
            Dict (i_α, i_β) → A tensor of shape (batch, total_dim, total_dim).
        """
        x_a2 = EvenSphere.patch_change_coords(x_alpha)
        x_b2 = EvenSphere.patch_change_coords(x_beta)
        patch_alpha_coords = {0: x_alpha, 1: x_a2}
        patch_beta_coords  = {0: x_beta,  1: x_b2}

        A_dict: Dict[Tuple[int, int], torch.Tensor] = {}
        for ia, ib in _PATCH_KEYS:
            x_combined = torch.cat([patch_alpha_coords[ia], patch_beta_coords[ib]], dim=1)
            A_dict[(ia, ib)] = self.patches[f'patch_{ia}{ib}'].forward_raw(x_combined)
        return A_dict


if __name__ == "__main__":
    print("Testing two-patch almost complex structure network...")

    dim = 4
    batch = 16
    sphere = EvenSphere(dimension=dim)

    net = TwoPatchAlmostComplexNet(
        manifold_dim=dim,
        hidden_layers=[64, 128, 64],
        activation='tanh',
        use_fourier_features=True,
    )

    x1 = sphere.sample_points_patch1(batch)
    J1, J2, x1, x2 = net(x1)
    print(f"J1 shape: {J1.shape},  J2 shape: {J2.shape}")

    I = torch.eye(dim).unsqueeze(0).expand(batch, -1, -1)
    err1 = torch.norm(torch.bmm(J1, J1) + I, dim=(1, 2)).mean()
    err2 = torch.norm(torch.bmm(J2, J2) + I, dim=(1, 2)).mean()
    print(f"||J1²+I||: {err1:.2e},  ||J2²+I||: {err2:.2e}")
    print("Test passed!")
