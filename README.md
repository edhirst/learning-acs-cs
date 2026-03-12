# Nijenhuis Tensor Minimisation via PINNs

Physics-informed neural network that learns an almost complex structure $J: TM \to TM$ minimising $\|N_J\|$. Supports two manifold types: `sphere` ($S^{2n}$) and `sphere_product` ($S^{m_\alpha} \times S^{m_\beta}$).

## Mathematics

An **almost complex structure** is a smooth bundle endomorphism satisfying $J^2 = -I$. Its integrability is obstructed by the **Nijenhuis tensor**:

$$N_J(X,Y) = [JX,JY] - J[JX,Y] - J[X,JY] + J^2[X,Y]$$

$N_J = 0$ everywhere $\iff$ $J$ is integrable (Newlander–Nirenberg). Among spheres, only $S^2$ and $S^6$ admit almost complex structures; whether $S^6$ admits an integrable one is an open problem.

## Architecture

Each patch sub-network (`PatchSubNet` in `src/models/architecture.py`):
1. **Fourier features**: $x \mapsto [\sin(2\pi B x),\, \cos(2\pi B x)]$
2. **MLP**: outputs $2n^2$ scalars parameterising $A(x)$ in the anti-commuting subspace $\{A \in \mathfrak{gl}(2n,\mathbb{R}) : AJ_0 + J_0 A = 0\}$
3. **Parametrisation**: $J(x) = e^{A(x)}\,J_0\,e^{-A(x)}$, so $J^2 = -I$ holds exactly for any weights

**`sphere`** — `TwoPatchAlmostComplexNet`: two `PatchSubNet`s covering $S^{2n}$ via stereographic ball patches.

**`sphere_product`** — `FourPatchProductNet`: four `PatchSubNet`s covering $S^{m_\alpha} \times S^{m_\beta}$ via patches $(i_\alpha, i_\beta) \in \{0,1\}^2$, each receiving the full concatenated coordinates. The per-patch $J_0$ is block-diagonal $J_0^\alpha \oplus J_0^\beta$ (with sign flips under patch transitions) when both factor dimensions are even.

## Loss Functions

The Nijenhuis tensor is assembled in chart coordinates via `torch.autograd.grad`:

$$N^k_{\;ij} = J^l_{\;i}\,\partial_l J^k_{\;j} - J^l_{\;j}\,\partial_l J^k_{\;i} + J^k_{\;l}(\partial_i J^l_{\;j} - \partial_j J^l_{\;i})$$

Overlap losses enforce cross-chart consistency $F(x) \cdot J \cdot F(x)^{-1} = J'$ at patch transitions, weighted to the equatorial region. All losses use radial Gaussian weighting with a flat plateau near the patch centre.

| Manifold | Nijenhuis terms | Overlap terms |
|---|---|---|
| `sphere` | 2 (one per patch) | 1 |
| `sphere_product` | 4 (one per patch pair) | 6 ($\binom{4}{2}$) |

## Setup

See [environment/](environment/) for setup instructions.

## Usage

```bash
# Train
python train.py --config configs/config_sphere2.yaml
python train.py --config configs/config_sphere_product.yaml

# Evaluate checkpoint
python evaluate.py outputs/<timestamp>/checkpoints/best_model.pt
```

## BibTeX Citation
```
raise NotImplementedError("Paper yet to be posted.")
```