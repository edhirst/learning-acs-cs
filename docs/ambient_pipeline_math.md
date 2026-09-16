# Unsupervised construction of almost complex and complex structures on sphere products: the mathematics of the ambient pipeline

This note documents the mathematics behind `extra/campaign_unit.py` (the unit that
`pbs/campaign_ambient.pbs` runs for every case). It covers the chart-free
ambient description of the geometry, the neural parametrizations, the losses,
the initialization and curriculum, the certification step, and the evaluation
protocol. Nothing in training uses a known solution; the classical structures
appear only in offline sanity checks, and those are flagged as such.

Cases: $S^2$, $S^4$, $S^6$, $S^1\times S^1$, $S^2\times S^2$, $S^3\times S^1$,
$S^3\times S^3$, $S^4\times S^2$, $S^5\times S^1$. Tasks: an almost complex
structure (ACS) and a complex structure (CS) on each.

---

## 1. Objects and goal

Let $M$ be a closed manifold of even dimension $d$. An **almost complex
structure** is a smooth endomorphism field $J$ of $TM$ with

$$J^2 = -\mathrm{Id}_{TM}.$$

Its **Nijenhuis tensor** is the $(1,2)$-tensor

$$N_J(U,V) = [JU,JV] - J[JU,V] - J[U,JV] + J^2[U,V] = [JU,JV] - J[JU,V] - J[U,JV] - [U,V],$$

and by Newlander–Nirenberg $J$ is induced by a **complex structure** (holomorphic
charts) if and only if $N_J \equiv 0$. Existence facts used to interpret the
results (never inside the optimisation):

* $S^{2n}$ admits an ACS iff $n = 1, 3$ (Borel–Serre); $S^4$ admits none.
* A product of even spheres is almost complex iff it is a product of copies of
  $S^2$, $S^6$, $S^2\times S^4$ (Sutherland); $S^1\times S^1$, $S^3\times S^1$,
  $S^3\times S^3$, $S^5\times S^1$ are parallelizable and carry complex
  structures (torus, Hopf, Calabi–Eckmann, Samelson).
* $S^6$ carries a complex structure (Alpöge–Claude / Engel, 2026), and hence so
  does $S^2\times S^4$ (Liang–Liu–Zhang, 2026, by blow-up, $\mu_2$-quotient and
  Atiyah flop; that structure has $c_1 = 0$).

The numerical goal is a **physics-informed neural network (PINN)** in the
strict sense: no data, no target field. The only training signals are the
residuals of the two defining equations, $J^2 + \mathrm{Id} = 0$ and $N_J = 0$,
sampled at random points of $M$. The pipeline must find (and certify) an ACS
on the eight cases where one exists, find a CS where one exists, and fail on
$S^4$.

---

## 2. The ambient-projector description (no charts, no frames)

Every case is a product of unit spheres

$$M = S^{m_1}\times\cdots\times S^{m_k} \subset \mathbb R^{m_1+1}\times\cdots\times\mathbb R^{m_k+1} = \mathbb R^N,\qquad d = \sum m_i,\quad N = \sum (m_i+1).$$

A point is $p = (x_1,\dots,x_k)$ with $|x_i| = 1$. The tangent space is
$T_pM = \bigoplus_i x_i^\perp$ and the orthogonal projector onto it is the
block-diagonal matrix

$$P(p) = \operatorname{diag}\big(I - x_1x_1^{\!\top},\ \dots,\ I - x_kx_k^{\!\top}\big),\qquad P^2 = P,\ P^{\top} = P,\ \operatorname{rank}P = d.$$

The induced metric is the round product metric, and $P$ is its identity
endomorphism on $TM$.

**Endomorphism fields.** A tangential endomorphism field is an $N\times N$
matrix field with $J = PJP$. Any matrix field $\tilde M(p)$ produces one by
$J = P\tilde M P$. This is how the network output is turned into a field of
endomorphisms of $TM$: the network never has to know about tangency.

**Tangent vector fields and brackets.** For $u\in\mathbb R^N$ the field
$U = Pu$ is tangent. For smooth ambient vector fields $U, V:\mathbb R^N\to\mathbb R^N$

$$[U,V] = (DV)\,U - (DU)\,V,\qquad (DV)_{ab} = \partial_b V_a,$$

and if $U,V$ are tangent along $M$ the result on $M$ is independent of the
ambient extension, so it is the Lie bracket on $M$. For a column field
$A(p)e_i$ of a matrix field $A$ one has $(D(Ae_i))\,w = \sum_k w_k\,\partial_kA\,e_i$;
everything is therefore expressible through $J$, $P$ and their ambient
Jacobians $DJ$, $DP$, all obtained by forward-mode automatic differentiation.

**Nijenhuis tensor on the overcomplete tangent frame.** Evaluate $N_J$ on the
$N$ tangent fields $U_i = Pe_i$:

$$N_{ij} := N_J(Pe_i, Pe_j) = [Je_i, Je_j] - J[Je_i, Pe_j] - J[Pe_i, Je_j] + J^2[Pe_i, Pe_j]\in\mathbb R^N .$$

(The last term is written with $J^2$, not $-P$, so that the same expression is
valid for the free two-term control on $S^4$, where $J^2\neq -P$.) Because
$\sum_i (Pe_i)(Pe_i)^{\!\top} = P$, the vectors $Pe_i$ form a Parseval (tight)
frame of $T_pM$, and for any bilinear $B$ on $T_pM$

$$\sum_{i,j=1}^N |B(Pe_i,Pe_j)|^2 = \sum_{a,b=1}^d |B(\varepsilon_a,\varepsilon_b)|^2$$

for any orthonormal basis $\varepsilon_a$ of $T_pM$. Hence

$$\|N_J\|^2(p) := \sum_{i,j}|N_{ij}(p)|^2$$

is exactly the intrinsic round-metric norm of the Nijenhuis tensor: coordinate
free, chart free, and comparable across cases. Offline check: for the
octonionic $G_2$-invariant ACS $J_{\rm oct}(x)v = x\times v$ on $S^6$ the code
returns $\|N\| = 19.5959$ (i.e. $\|N\|^2 = 384$) at every point, constant to
$10^{-16}$, as the homogeneity of that structure demands; the octonionic-product
ACS on $S^4\times S^2$ gives $15.17$.

**Why this replaces charts.** The original chart formulation carried two costs:
a cross-chart gluing loss that competes with the Nijenhuis loss, and chart
coordinates that inflate derivatives ($\|\partial J\|\approx 34$ for the exact
Calabi–Eckmann structures, which is what defeats gradient training of a
first-order residual). The ambient field is one global object, so there is
nothing to glue, and derivatives are those of the round geometry.

**Sampling.** Points are sampled from the round product measure: each factor
$x_i = g/|g|$ with $g$ standard Gaussian in $\mathbb R^{m_i+1}$. All losses are
Monte-Carlo expectations over this measure with mini-batches of 192 points.

**Retraction and tangential gradient** (used by the adversarial step): for a
point $p$ and an ambient vector $g$,
$\operatorname{grad}_T g = (g_i - \langle g_i, x_i\rangle x_i)_i$ and
$\operatorname{retract}(p + \delta) = \big((x_i+\delta_i)/|x_i+\delta_i|\big)_i$.

*(A global-frame formulation, valid only on the four parallelizable cases, is
retained in the code as `--formulation frame` but is not used by the campaign.)*

---

## 3. Neural parametrizations

All networks take the ambient point $p\in\mathbb R^N$ as input and return an
$N\times N$ matrix (or several). Two families are used.

**Polynomial fields (harmonic classes).** $\tilde M(p) = \sum_\alpha W_\alpha\,\phi_\alpha(p)$ with $W_\alpha\in\mathbb R^{N\times N}$ and features from one of the classes

| class | features $\phi_\alpha(p)$ | harmonic content on $M$ |
|---|---|---|
| `c` | $1$ | degree 0 |
| `lin` | $p_a$ | degree 1 |
| `bilin` | $x_{A,a}\,x_{B,b}$ (two-factor products only) | bidegree $(1,1)$ |
| `quad_h` | $p_ap_b$, $a\le b$ | degrees 2 and 0 |
| `lin+c`, `lin+bilin`, `all2` | unions of the above; `all2` $=\{1, p_a, p_ap_b\}$ | degree $\le 2$ |

The restriction of a homogeneous polynomial of degree $\ell$ on $\mathbb R^{m+1}$
to $S^m$ is a sum of spherical harmonics of degrees $\ell, \ell-2, \dots$, so the
ladder is literally "lowest harmonics first". Weights are initialised
$W_\alpha\sim\mathcal N(0, \sigma^2/n_{\rm feat})$ with $\sigma = 0.3$.

**MLP.** $\tilde M(p) = \mathrm{reshape}(\mathrm{MLP}_\theta(p))$, GELU
activations, widths $256^{\times 4}$ for the six-dimensional cases and
$128^{\times 3}$ otherwise; the last layer is scaled by $\sigma$ ($0.3$ for
ACS-refinement, $0.5$ for the CS nets, $0.02$ when the MLP is an additive
refinement of a polynomial field).

The composite ACS field is $\tilde M = \mathrm{poly}(p) + \mathrm{mlp}(p)$
(`CurriculumNet`), the CS nets are plain MLPs.

---

## 4. The ACS stage

### 4.1 Free field and loss

The candidate is the free tangential field

$$K(p) = P(p)\,\tilde M_\theta(p)\,P(p),\qquad r(p) = \|K^2 + P\|_F,\qquad \mathcal L_{\rm ACS}(\theta) = \mathbb E_p\, r(p)^2 .$$

$r = 0$ at a point iff $K|_{T_pM}$ is an almost complex structure there. No base
point, no reference structure.

### 4.2 Why the plain flow fails: orientation domains and defects

At a point, the set of ACS on $T_pM$ has two connected components (the two
orientations), each homotopy equivalent to $SO(d)/U(d/2)$; for $d = 6$ this is
$\mathbb{CP}^3$, with $\pi_2 = \mathbb Z$. A randomly initialised $K$ chooses
its orientation locally, so gradient descent produces domains of $+J$ and $-J$
separated by walls on which $K$ degenerates (and, in dimension six, additional
codimension-three vortex-type defects). These configurations are **strict local
minima** of the mean loss, at every resolution. Explicitly on $S^2$, for the
two-parameter family $K = (a + b\,f)J_{\rm std}$ with $f = \langle x, n\rangle$
(uniform on $[-1,1]$ under the round measure)

$$R(a,b) = 2\,\mathbb E\big[(1-(a+bf)^2)^2\big],\qquad R(0,b) = 2\Big(1 - \tfrac{2b^2}{3} + \tfrac{b^4}{5}\Big),$$

minimised at $b^2 = 5/3$, where $\partial_a^2R = 2(-4 + 4b^2) > 0$: the wall
configuration (one hemisphere $+J$, the other $-J$) is a local minimum in the
orientation-mixing direction $a$. This is what the earlier campaign observed on
every non-parallelizable case (residual sup $\approx\sqrt 2 = \|P\|_F$ at the
walls, i.e. $K\approx 0$ there).

### 4.3 Harmonic ladder (the curriculum)

The remedy is a curriculum in the *function class*, not in the loss: train the
polynomial field in the lowest harmonic class first and climb

$$\texttt{c}\ \to\ \texttt{lin}\ \to\ \texttt{bilin}\ \to\ \texttt{quad\_h}\ \to\ \texttt{lin+c}\ \to\ \texttt{lin+bilin}\ \to\ \texttt{all2},$$

accepting the first class whose trained field passes certification (§4.5). A
field of low harmonic degree cannot contain a small defect core, so the flow in
such a class either settles on one global structure or stays loudly far from
$r = 0$. Each class is trained for 3000 Adam steps (lr $5\cdot10^{-3}$, cosine
decay) with 192-point batches and costs seconds to a minute.

Where an ACS exists it is typically found in the lowest classes: the accepted
classes in validation were `c` on the frame cases, `lin` on $S^2$, $S^2\times S^2$,
$S^6$, $S^3\times S^1$, `bilin` on $S^4\times S^2$, $S^3\times S^3$, $S^5\times S^1$,
$S^1\times S^1$. That the classical structures are themselves of low degree
($x\times v$ on $S^2$ and $S^6$ is linear, the octonionic-product ACS on
$S^4\times S^2$ is bilinear) is an a-posteriori explanation, not an input: the
ladder is a generic "low frequencies first" device, and on $S^6$ it also finds
non-octonionic exact ACS from other seeds.

### 4.4 Dirichlet homotopy (retry only)

If a class is rejected, it is retried once with an annealed smoothness term,
$\mathcal L = \mathcal L_{\rm ACS} + \lambda(t)\,\mathcal E_D$, with the
**intrinsic** Dirichlet energy

$$\mathcal E_D(K) = \mathbb E_p \sum_{a=1}^N \big\|P\,(D_{Pe_a}K)\,P\big\|_F^2 = \mathbb E_p\,\|\nabla K\|^2 ,$$

where $D_{Pe_a}K = \sum_k P_{ka}\partial_kK$ and $\nabla$ is the Levi-Civita
connection of the round metric (for a submanifold, $\nabla_XJ = P\,(D_XJ)\,P$ for
tangential $J$; the tight-frame identity makes the sum over ambient directions
intrinsic). $\lambda(t)$ decreases geometrically from $3$ to $10^{-3}$ over the
class's training steps.

Mechanism: on the cases whose structure is parallel ($\nabla J = 0$: $S^2$,
$S^2\times S^2$, the flat torus) a wall carries derivative energy while the
structure does not, so a large $\lambda$ shrinks the wall amplitude $b$ until
the orientation direction $a$ becomes unstable ($\partial_a^2R<0$ once $b^2<1$)
and the flow picks one global orientation; the anneal then removes the term.
On $S^6$ and $S^4\times S^2$ no ACS is parallel ($\mathcal E_D(J_{\rm oct}) =
24$, while the wall configuration $\langle x,n\rangle J_{\rm oct}$ has $60/7$),
so the term must not be applied there; the ladder never reaches the retry on
those cases because the plain attempt already certifies. A scale-invariant
variant (energy of the orthogonal polar factor of $K$) was tried and rejected:
it collapses $K\to 0$.

### 4.5 Refinement and certification

After acceptance the MLP is switched on as an additive correction and both
parts are trained on $\mathcal L_{\rm ACS}$ (Adam lr $10^{-3}$, cosine to
$10^{-5}$, gradient clipping at 1, 10 000 steps), followed by an L-BFGS polish.

**Lemma (degeneracy bound).** If $K|_{T_pM}$ has a real eigenvalue
$\lambda$ with unit eigenvector $v\in T_pM$, then $(K^2+P)v = (\lambda^2+1)v$,
so $\|K^2+P\|_{\rm op}\ge 1$ and $r(p)\ge 1$. Conversely, if $r(p) < 1$ then
$K|_{T_pM}$ has no real eigenvalue, hence is conjugate to an ACS.

**Adversarial sup.** Existence is therefore a *sup* statement, and a sampled
sup is not enough: on $S^4$ any continuous $K$ must degenerate somewhere, but
the degenerate set can be tiny (a trained $S^4$ field had sampled sup $0.35$ on
$6000$ points and a true sup of $2$). The certificate used is

$$\hat r_{\rm adv} = \max_{\rm starts}\ \text{Riemannian ascent of } r,$$

with starts at the 256 worst points of an $8192$-point sample and 150 ascent
steps $p\leftarrow\operatorname{retract}(p + \eta\,\operatorname{grad}_T r/|\operatorname{grad}_T r|)$
with a per-start line search ($\eta\times 1.2$ on increase, $\times 0.5$ otherwise).
The field is **admissible** iff $\hat r_{\rm adv} < 1$. This is much stronger
than sampling but still not a proof; the paper states it as evidence.

**Polar correction (exactness).** On an admissible field set
$S = -K^2 + (I - P)$. On $T_pM$ its spectrum is $\{-\lambda^2\}$ over the
(non-real) eigenvalues $\lambda = \alpha\pm i\beta$, $\beta\neq0$, i.e.
$\beta^2-\alpha^2 \mp 2i\alpha\beta$, never on the negative real axis, so the
principal inverse square root exists and

$$J^* = K\,S^{-1/2},\qquad (J^*)^2 = K^2S^{-1} = -P\ \text{ on } T_pM,\quad J^* = 0\ \text{ on } (T_pM)^\perp ,$$

because $K$ commutes with $S$. $J^*$ is an exact ACS ($\|J^{*2}+P\| \sim 10^{-10}$
in double precision) homotopic to $K$. $S^{-1/2}$ is computed by the coupled
Newton–Schulz iteration on $S/2$,

$$Y_0 = S/2,\ Z_0 = I,\qquad T_k = \tfrac12(3I - Z_kY_k),\ Y_{k+1} = Y_kT_k,\ Z_{k+1} = T_kZ_k,\qquad Z_k\to (S/2)^{-1/2},$$

which converges quadratically whenever the spectrum of $S$ lies in the disc
$|s-2|<2$, guaranteed by $r<1$ (then every eigenvalue of $S$ is within distance
$<1$ of $1$). Unlike an eigendecomposition it is built from matrix products
only, so it is differentiable, forward-AD and `vmap` safe; this is what allows
$J^*$ to serve as a base point in the CS stage. (The first version of the code
used a symmetric eigensolver on $S$, which is wrong because $S$ is not
symmetric for a non-normal $K$; the general eigendecomposition is kept as a
cross-check.)

**Verdict.** *Certified exact ACS* iff $\hat r_{\rm adv}(K) < 1$ and
$\hat r_{\rm adv}(J^*) < 10^{-8}$. On $S^4$ every class and every retry is
rejected ($\hat r_{\rm adv}\approx\sqrt 2$ – $2$), which is the topological
obstruction showing up as a loud, resolution-independent failure of the sup,
even though the *mean* residual can be made small by concentrating the defect.

---

## 5. The CS stage

### 5.1 Exact-ACS parametrization

With $J_{\rm base} = J^*$ the unit's own certified ACS (frozen), the candidate is

$$J = e^{B_2}\,e^{B_1}\,J_{\rm base}\,e^{-B_1}\,e^{-B_2},\qquad B_s = P\,\mathrm{mlp}_s(p)\,P .$$

Since $B_s$ commutes with $P$, so does $e^{B_s}$, and $J^2 = e^{B}J_{\rm base}^2e^{-B} = -P$
exactly: the algebraic constraint is built in and the only equation left is
$N_J = 0$. Two conjugations suffice to reach every ACS field of the same
orientation: two ACS at a point are conjugate by some $g\in GL^+(d)$, and the
polar decomposition $g = e^{\Sigma}R$ with $\Sigma$ symmetric and $R\in SO(d) =
\exp(\mathfrak{so}(d))$ writes $g$ as a product of two exponentials (a single
exponential is not onto $GL^+$). The $B_s$ are *general* matrices: restricting
them to be antisymmetric (orthogonal conjugation) would confine the search to
round-orthogonal structures, which LeBrun's theorem excludes on $S^6$.

Because $N_{-J} = N_J$, the orientation of the base point is immaterial.

### 5.2 Initialization and loss

The CS nets are initialised at $O(1)$ ($\sigma = 0.5$), so the search starts
from a *generic non-integrable* structure: initial $\mathbb E\|N_J\|^2$ is
$10^3$–$10^6$ in the campaign. The loss is the mean squared intrinsic
Nijenhuis norm of §2,

$$\mathcal L_{\rm CS}(\theta) = \mathbb E_p \sum_{i,j}|N_{ij}(p)|^2 ,$$

with $DJ$ by forward-mode AD through the exponentials. Optimisation: Adam
(lr $10^{-3}$, cosine to $10^{-5}$, clipping at 1, batch 192; 20k / 40k / 80k
steps for $d = 2 / 4 / 6$), optional residual-adaptive resampling (half of each
batch from the worst of $2\times$ candidates), then L-BFGS (strong-Wolfe line
search) on a fixed 768-point batch.

No certified base means no CS unit: on $S^4$ the stage refuses, so the
pipeline's failure on $S^4$ is inherited from a topological fact, which is the
only kind of numerical failure that carries information (a small Nijenhuis
residual is evidence of existence, but a large one is never evidence of
non-existence).

In real dimension two $N_J\equiv 0$ for every $J$, so the CS stage on $S^2$ and
$S^1\times S^1$ is trivially exact.

---

## 6. Evaluation

Each unit reports, on fresh points:

* dense mean and max of $r$ and of $\|N_J\|$;
* the adversarial sups $\hat r_{\rm adv}$ and $\widehat{\|N\|}_{\rm adv}$ (same
  ascent applied to the Nijenhuis norm);
* for ACS units, the corrected field's residuals and the verdict; the accepted
  ladder class and the full ladder log;
* the paper's original **chart protocol**, for comparability with the earlier
  table: the global field is pulled back to the ball charts
  $\sigma$ of each patch by $J_{\rm chart} = (D\sigma)^{+}\,J\,D\sigma$
  (Moore–Penrose pseudoinverse), then the hard-cutoff per-patch Nijenhuis mean
  (CS column) and the cross-patch overlap discrepancy (ACS column) are computed
  exactly as before. For a single global field the overlap column is
  identically $\sim 10^{-10}$, which is why the ACS quality is read from
  $\hat r_{\rm adv}$ and the certificate instead.

**Chern class check (offline).** For a tangential ACS $J$ the $(1,0)$-bundle is
the image of the complex projector $\Pi = (P + iJ)/2\subset\mathbb C^N$, whose
projected connection has curvature $\Pi\,d\Pi\wedge d\Pi$, so on a 2-cycle
$\Sigma$

$$\langle c_1(T^{1,0}M), [\Sigma]\rangle = \frac{i}{2\pi}\int_\Sigma \operatorname{tr}\big(\Pi\,[\partial_\theta\Pi,\partial_\varphi\Pi]\big)\,d\theta\,d\varphi .$$

$c_1$ is a homotopy invariant of the ACS and is constant along the conjugation
flow of §5, so the CS search can only reach complex structures in the base
point's Chern class. The check (validated on $S^2$ and $S^2\times S^2$, $c_1 = \pm 2$)
gives $c_1 = 0$ for the $S^4\times S^2$ base point, the class of the
Liang–Liu–Zhang structure: no topological obstruction remains on that row.

---

## 7. The algorithm in one page

```
for each case M, seed:
  # ---- ACS unit -------------------------------------------------------
  for class in [c, lin, bilin, quad_h, lin+c, lin+bilin, all2]:
      for attempt in [plain, plain + annealed Dirichlet (only if plain rejected)]:
          K = P poly_class(p) P ; minimise E||K^2 + P||^2  (3000 Adam steps)
          if adversarial_sup(||K^2+P||) < 1: accept; break
  K = P (poly + mlp)(p) P ; minimise E||K^2 + P||^2  (10k Adam + L-BFGS)
  if adversarial_sup(||K^2+P||) < 1:
      J* = K (-K^2 + I - P)^{-1/2}   (Newton–Schulz)      -> certified exact ACS
  else: NOT certified (S^4)
  # ---- CS unit ----------------------------------------------------------
  if no certified ACS: refuse
  J = e^{B2} e^{B1} J* e^{-B1} e^{-B2},  B_s = P mlp_s(p) P, O(1) init
  minimise E ||N_J||^2  (Adam + L-BFGS) ; report dense + adversarial residuals
  # ---- scoring ---------------------------------------------------------
  chart-protocol numbers (for the table), c1 check offline
```

---

## 8. Numbers worth remembering (offline checks, not inputs)

| quantity | value |
|---|---|
| $\|N\|$ of the octonionic ACS on round $S^6$ | $19.5959$ (constant; $\|N\|^2 = 384$) |
| $\|N\|$ of the octonionic-product ACS on $S^4\times S^2$ | $15.17$ |
| intrinsic Dirichlet energy $\mathcal E_D$: $S^2$, $S^2\times S^2$, torus, frame $J_0$ | $0$ (parallel) |
| $\mathcal E_D$ of the octonionic ACS on $S^6$ and of the $S^4\times S^2$ ACS | $24$ |
| wall configuration $\langle x,n\rangle J_{\rm oct}$ on $S^6$: $\mathcal E_D$ | $60/7\approx 8.6$ (why the homotopy fails there) |
| $c_1[S^2]$ of the $S^4\times S^2$ base ACS | $0$ |
| residual at a degenerate point | $\ge 1$; at $K = 0$: $\|P\|_F = \sqrt d$ |

---

## 9. Files

| file | content |
|---|---|
| `src/geometry/ambient.py` | `AmbientManifold`: sampling, projector $P$ |
| `src/loss/ambient_loss.py` | ambient brackets, $N_{ij}$, $\|N\|$, $DJ$/$DP$ by forward AD |
| `src/geometry/explicit_acs.py` | case registry; classical ACS (offline checks / optional bases) |
| `extra/campaign_unit.py` | the unit: ladder, Dirichlet retry, MLP refinement, adversarial sup, Newton–Schulz polar, CS stage, chart protocol |
| `extra/campaign_summary.py` | before/after table with verdicts |
| `extra/chern_class_ambient.py` | $c_1$ on the $S^2$ cycle |
| `extra/acs_homotopy_probe.py` | the experiments behind §4.2–4.4 |
| `pbs/campaign_ambient.pbs`, `jobs_ambient_{acs,cs}.txt` | the HPC job (phase A: ACS, phase B: CS from own ACS) |

## References

Newlander–Nirenberg (1957); Borel–Serre (1953); Sutherland (1965); LeBrun,
*Orthogonal complex structures on $S^6$* (1987); Calabi–Eckmann (1953);
Samelson (1953); Alpöge–Claude / Engel, *A compact complex threefold fibred by
tori over the projective line, and the six-sphere* (2026); Liang–Liu–Zhang,
*Remarks on the complex structures on $\mathbb P^3$ and $S^2\times S^4$*
(arXiv:2608.29651); Higham, *Functions of Matrices* (2008) for Newton–Schulz;
Wall (1966) for the classification used by Liang–Liu–Zhang.
