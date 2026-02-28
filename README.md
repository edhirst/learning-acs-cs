# Nijenhuis Tensor Minimisation via PINNs

Physics-informed neural network that learns an almost complex structure $J: TM \to TM$ on even-dimensional spheres $S^{2n}$ minimising $\|N_J\|$.

## Mathematics

An **almost complex structure** is a smooth bundle endomorphism satisfying $J^2 = -I$. Its integrability is obstructed by the **Nijenhuis tensor**:

$$N_J(X,Y) = [JX,JY] - J[JX,Y] - J[X,JY] + J^2[X,Y]$$

$N_J = 0$ everywhere $\iff$ $J$ is integrable (Newlander–Nirenberg). Among spheres, only $S^2$ and $S^6$ admit almost complex structures; whether $S^6$ admits an integrable one is an open problem.

## Architecture

`TwoPatchAlmostComplexNet` wraps two independent `PatchSubNet` instances. Given patch-1 coordinates $x_1 \in B(0,1)$, it computes $x_2 = \varphi(x_1)$ and evaluates both sub-networks:

```
x1 ──► PatchSubNet1 ──► J1
x1 ──► φ ──► x2 ──► PatchSubNet2 ──► J2
```

Each `PatchSubNet`:
1. **Fourier features**: $x \mapsto [\sin(2\pi B x),\, \cos(2\pi B x)]$
2. **MLP**: outputs $\dim^2$ scalars forming an unconstrained matrix $A(x) \in \mathfrak{gl}(2n,\mathbb{R})$
3. **Cayley-like parametrisation**: $J(x) = e^{A(x)}\,J_0\,e^{-A(x)}$

$J^2 = -I$ holds exactly for any network weights since $(e^A J_0 e^{-A})^2 = e^A J_0^2 e^{-A} = -I$.

## Loss Functions

**Per-patch Nijenhuis** (`NijenhuisTensorLoss`): assembles the coordinate expression

$$N^k_{\;ij} = J^l_{\;i}\,\partial_l J^k_{\;j} - J^l_{\;j}\,\partial_l J^k_{\;i} + J^k_{\;l}(\partial_i J^l_{\;j} - \partial_j J^l_{\;i})$$

using $\partial J/\partial x$ obtained by `torch.autograd.grad` in the ball chart.

**Overlap loss** (`OverlapLoss`): enforces cross-chart consistency

$$L_\text{overlap} = \|F(x_1)J_1 F(x_1)^{-1} - J_2\|^2 + \|F(x_2)J_2 F(x_2)^{-1} - J_1\|^2$$

**Total loss**: $L = \lambda_\text{nij}(\|N_{J_1}\| + \|N_{J_2}\|) + \lambda_\text{overlap}\,L_\text{overlap}$

## Setup

```bash
conda env create -f environment/environment.yml
conda activate nijenhuis
```

Or with pip:

```bash
pip install -r environment/requirements.txt
```

## Usage

```bash
# Train on sphere (default config = S²)
python train.py --config configs/config_sphere2.yaml

# Evaluate checkpoint
python evaluate.py outputs/<timestamp>/checkpoints/best_model.pt
```

## BibTeX Citation  
``` 
raise NotImplementedError("Paper yet to be posted.")
```