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
    All patch combinations share the same full-dimensional canonical structure
    J₀(d), d = m_α + m_β, regardless of whether the individual factor dimensions
    are even or odd.

In all cases the exp(A) J₀ exp(−A) parametrisation guarantees J² = −I exactly
for any network weights.
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Dict, List, Optional, Tuple

from src.geometry.sphere import EvenSphere


def _taylor_matrix_exp(A: torch.Tensor, num_terms: int) -> torch.Tensor:
    """
    Matrix exponential via truncated Taylor series with scaling-and-squaring.

    Computes (∑_{k=0}^{num_terms} B^k/k!)^{2^s} where B = A/2^s and
    s = max(0, ⌈log₂ ‖A‖₁⌉) ensures ‖B‖₁ ≤ 1.

    Args:
        A: Input matrices, shape (batch, d, d).
        num_terms: Number of Taylor terms K (not counting the k=0 identity).

    Returns:
        Approximation of exp(A), shape (batch, d, d).
    """
    if A.device.type == 'mps':
        return _taylor_matrix_exp(A.cpu(), num_terms).to(A.device)

    batch, d, _ = A.shape
    dtype, device = A.dtype, A.device

    # Scale so ‖B‖₁ ≤ 1 to improve series convergence.
    norm = torch.linalg.matrix_norm(A, ord=1).amax()
    s = int(torch.clamp(torch.ceil(torch.log2(norm.clamp(min=1.0))), min=0).item())
    B = A / (2.0 ** s) if s > 0 else A

    I = torch.eye(d, dtype=dtype, device=device).unsqueeze(0)
    result = I.expand(batch, -1, -1).clone()
    term   = result.clone()
    for k in range(1, num_terms + 1):
        term   = torch.bmm(term, B) / k
        result = result + term

    for _ in range(s):
        result = torch.bmm(result, result)

    return result


def _matrix_exp(A: torch.Tensor, num_terms: Optional[int] = None) -> torch.Tensor:
    """Matrix exponential; uses torch.matrix_exp (Padé-13) by default.

    Args:
        A: Input matrices, shape (batch, d, d).
        num_terms: If None, uses torch.matrix_exp.  If an integer K, uses a
            K-term Taylor series with scaling-and-squaring instead.
    """
    if num_terms is not None:
        return _taylor_matrix_exp(A, num_terms)
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
        B = torch.randn(input_dim, num_features, dtype=torch.float64) * scale
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


class ResidualBlock(nn.Module):
    """
    Residual block: two linear layers with an activation and a skip connection.

    The residual path is: act(linear2(act(linear1(x)))).
    The skip connection is the identity when in_dim == out_dim, or a linear
    projection (bias=False) otherwise.

    Args:
        in_dim:  Input feature dimension.
        out_dim: Output feature dimension.
        act_cls: Activation class (e.g. nn.GELU).
    """

    def __init__(self, in_dim: int, out_dim: int, act_cls: type):
        super().__init__()
        self.linear1 = nn.Linear(in_dim, out_dim)
        self.act1    = act_cls()
        self.linear2 = nn.Linear(out_dim, out_dim)
        self.act2    = act_cls()
        self.skip    = nn.Identity() if in_dim == out_dim else nn.Linear(in_dim, out_dim, bias=False)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.act2(self.linear2(self.act1(self.linear1(x))) + self.skip(x))


class ResidualMLP(nn.Module):
    """
    MLP built from residual blocks followed by a linear output layer.

    Has the same call interface as nn.Sequential so it can replace
    self.network without touching forward() or forward_raw().

    Args:
        input_dim:     Input feature dimension.
        hidden_layers: List of hidden layer widths; one ResidualBlock per entry.
        act_cls:       Activation class.
        output_dim:    Dimension of the final linear projection.
    """

    def __init__(
        self,
        input_dim: int,
        hidden_layers: List[int],
        act_cls: type,
        output_dim: int,
    ):
        super().__init__()
        blocks: List[nn.Module] = []
        prev = input_dim
        for h in hidden_layers:
            blocks.append(ResidualBlock(prev, h, act_cls))
            prev = h
        self.blocks = nn.ModuleList(blocks)
        self.output = nn.Linear(prev, output_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        for block in self.blocks:
            x = block(x)
        return self.output(x)


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
        matrix_exp_terms: Optional[int] = None,
        use_residual_blocks: bool = False,
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
            matrix_exp_terms: If None, uses torch.matrix_exp (Padé-13 approximant).
                If a positive integer K, uses a K-term Taylor series with
                scaling-and-squaring instead.
            use_residual_blocks: If True, builds the MLP from ResidualBlocks
                instead of plain linear layers.
        """
        super().__init__()

        if manifold_dim % 2 != 0:
            raise ValueError("manifold_dim must be even")

        self.manifold_dim = manifold_dim
        self.complex_dim  = manifold_dim // 2
        self.matrix_exp_terms = matrix_exp_terms
        # Network outputs all d² components of A ∈ gl(d,ℝ) directly.
        # J = exp(A) J₀ exp(−A) satisfies J² = −I for any A; no subspace
        # restriction is required.
        self.num_dof = manifold_dim * manifold_dim

        # Canonical complex structure: use the supplied matrix when given,
        # otherwise build the standard block-diagonal J₀ (with sign).
        if j0_matrix is not None:
            J0 = j0_matrix
        else:
            J0 = j0_sign * self._build_J0(manifold_dim)
        self.register_buffer('J0', J0)

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
        if use_residual_blocks:
            self.network = ResidualMLP(input_dim, hidden_layers, act_cls, self.num_dof)
        else:
            layers = []
            prev = input_dim
            for h in hidden_layers:
                layers += [nn.Linear(prev, h), act_cls()]
                prev = h
            layers.append(nn.Linear(prev, self.num_dof))
            self.network = nn.Sequential(*layers)

        self._initialize_weights(initialization)

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
        batch = params.shape[0]
        A = params.reshape(batch, self.manifold_dim, self.manifold_dim)

        expA  = _matrix_exp(A,  self.matrix_exp_terms)
        expnA = _matrix_exp(-A, self.matrix_exp_terms)
        J0 = self.J0.unsqueeze(0).expand(batch, -1, -1)
        return expA @ J0 @ expnA

    def forward_raw(self, x: torch.Tensor) -> torch.Tensor:
        """Return A(x) ∈ gl(d,ℝ): the generator matrix before the exp transform.

        Args:
            x: Ball-patch coordinates, shape (batch, manifold_dim).

        Returns:
            A: Shape (batch, manifold_dim, manifold_dim).
        """
        features = self.fourier(x) if self.fourier is not None else x
        params = self.network(features)                      # (batch, num_dof)
        return params.reshape(params.shape[0], self.manifold_dim, self.manifold_dim)


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
        matrix_exp_terms: Optional[int] = None,
        use_residual_blocks: bool = False,
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
            matrix_exp_terms: Passed to PatchSubNet; see its docstring.
            use_residual_blocks: Passed to PatchSubNet; see its docstring.
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
            matrix_exp_terms=matrix_exp_terms,
            use_residual_blocks=use_residual_blocks,
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

    # All patches use the full-dimensional canonical structure J₀(d).
    J0_full = PatchSubNet._build_J0(d)
    for ia, ib in _PATCH_KEYS:
        j0s[(ia, ib)] = J0_full

    # [archived] Block-diagonal initialisation for both-even factor dimensions:
    # if dim_alpha % 2 == 0 and dim_beta % 2 == 0:
    #     # j0^{(ia,ib)} = block( (-1)^{ia} J₀^α, (-1)^{ib} J₀^β )
    #     # Reflects det(F_eq) = -1 for even-dimensional factors.
    #     J0a = PatchSubNet._build_J0(dim_alpha)   # (m_α, m_α)
    #     J0b = PatchSubNet._build_J0(dim_beta)    # (m_β, m_β)
    #     for ia, ib in _PATCH_KEYS:
    #         J0 = torch.zeros(d, d)
    #         J0[:dim_alpha, :dim_alpha] = ((-1) ** ia) * J0a
    #         J0[dim_alpha:, dim_alpha:] = ((-1) ** ib) * J0b
    #         j0s[(ia, ib)] = J0

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
        matrix_exp_terms: Optional[int] = None,
        use_residual_blocks: bool = False,
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
            matrix_exp_terms: Passed to PatchSubNet; see its docstring.
            use_residual_blocks: Passed to PatchSubNet; see its docstring.
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
            matrix_exp_terms=matrix_exp_terms,
            use_residual_blocks=use_residual_blocks,
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
