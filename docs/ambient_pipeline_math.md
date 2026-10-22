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

The implementation is not limited to two factors. The additional registered
case $S^4\times S^1\times S^1$ exercises the same ambient pipeline with three
factors ($d=6$, $N=9$); it is kept separate from the original nine-case table.

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

Sampling, projection, ambient Nijenhuis evaluation, retraction, adversarial
search, and the final $2^k$-chart cross-check all iterate over this factor
list. Thus adding a case is registry/configuration work rather than a new
geometric implementation.

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
| `bilin` | $x_{A,a}\,x_{B,b}$ for every pair of distinct factors $A<B$ | pairwise multidegree $(1,1)$ |
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

$$\texttt{c}\ \to\ \texttt{lin}\ \to\ \texttt{bilin}\ \to\ \texttt{quad\_h}\ \to\ \texttt{lin+c}\ \to\ \texttt{lin+bilin}\ \to\ \texttt{all2}\ \to\ \texttt{all2+trilin},$$

(the last rung, degree $\le 2$ plus the three-factor products
$x_{A,a}x_{B,b}x_{C,c}$, exists only for products of $\ge 3$ factors; it is the
multidegree-$(1,1,1)$ harmonic content and carries the frame ACS of
$S^4\times S^1\times S^1$, see `extra/test_poly_cs.py`)

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

### 4.6 Orientation control (third attempt)

If a class is rejected by the plain and the Dirichlet attempts, it is retried
with a **consistent orientation penalty**. Let $E(p)$ be an orthonormal frame of
$T_pM$, positively oriented for the global orientation in which
$(n_1,\dots,n_k,E)$ is a positive basis of $\mathbb R^N$ ($n_A$ the outward
normal of factor $A$), and

$$o(p) = \frac{\operatorname{Pf}\big(E^{\top}(K^{\top}-K)E\big)}{2^{d/2}} .$$

**Lemma.** If $A$ is an ACS on $\mathbb R^d$ then $\omega = A^{\top} - A$ is
nondegenerate, $|\operatorname{Pf}\omega| \ge 2^{d/2}$, and
$\operatorname{sign}\operatorname{Pf}\omega$ is the complex orientation of $A$.
*Proof.* $h = I + A^{\top}A$ is positive and $A$-invariant
($h(Au,Av) = h(u,v)$ because $A^2 = -I$), and $\omega(u,v) = h(Au, v)$, so
$\omega$ is an $A$-compatible symplectic form ($\omega(u,Au) = h(Au,Au) > 0$) and
$\det\omega = \det h = \prod(1+\sigma_i^2)$; the singular values of $A$ come in
pairs $\sigma, \sigma^{-1}$ ($A^{-1} = -A$), so $\det h \ge 4^{d/2}$. $\square$

The penalty $\mu\,\mathbb E\,\mathrm{relu}(1 - s\,o)^2$, with $s$ the sign of
$\mathbb E\,o$ at initialisation (the majority orientation), therefore vanishes
on *every* ACS of that orientation and is $\approx\mu$ on the $K\approx 0$ cores of
orientation walls (§4.2): it removes the wall minima without selecting a
particular structure. Fixing the orientation loses nothing on these manifolds
(every product of spheres has an orientation-reversing isometry). Validated:
$|o| \ge 1$ with one sign everywhere on every certified ACS, and on
$S^4\times S^1\times S^1$ it is what certifies (plain and Dirichlet attempts stay
at the wall value $\hat r_{\rm adv}\approx 1.42$; the orientation attempt gives
$\hat r_{\rm adv} = 0.12$).

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

### 5.3 Polynomial free-field search with a Levenberg–Marquardt finisher (`--cs-mode poly`)

The conjugation parametrization of §5.1 has to express an integrable structure
as $gJ_{\rm base}g^{-1}$ with $J_{\rm base}$ an MLP-refined learned field, so $g$
is complicated even when the integrable structure is simple. The poly mode
searches directly in a small class instead:

$$K = P\,M_W(p)\,P,\qquad M_W(p) = \sum_f \phi_f(p)\,W_f,$$

with $\phi$ all monomials of total degree $\le k$ (optionally plus the
three-factor products, class `2t`), the same "lowest harmonics first" idea as
the ACS ladder. $K$ is linear in $W$, so $N_K$ is a *quadratic* map of $W$ and
the problem in a fixed class is a small nonlinear least-squares problem,

$$\min_W\ \mathbb E\,\|N_K\|^2 + w\,\mathbb E\,\|K^2+P\|^2\qquad (w = 10),$$

zero exactly at an integrable ACS of the class. Initialisation is the
least-squares projection of the unit's own certified $J^*$ into the class (it
must itself be admissible); the result is certified and polar-corrected by §4.5,
and $N$ is reported for the corrected exact ACS $J^*$. **Offline check (not an
input):** the Hopf structures on $S^3\times S^1$, $S^5\times S^1$
($J\partial_t = \xi$, $J\xi = -\partial_t$, $J|_H = i$) and the Calabi–Eckmann
structure on $S^3\times S^3$ ($J\xi_1 = \xi_2$, $J|_{H_1\oplus H_2} = i$) are
exactly $P M P$ with $M$ of degree 2 (fit residual $10^{-15}$, $N = 0$), and the
product structure on $S^2\times S^2$ is linear. So on these cases the class
contains an exact solution, and failure would be an optimisation failure.

*Residual layout.* With $E$ an orthonormal tangent frame, $N$ is
tangent-valued and antisymmetric, so
$\sum_{i,j}|N(Pe_i,Pe_j)|^2 = 2\sum_{a<b}|E^{\top}N(E_a,E_b)|^2$ and
$K^2 + P = E(E^{\top}K^2E + I_d)E^{\top}$: $d\binom d2 + d^2$ rows per point
(126 for $d=6$) instead of $N^3+N^2$, with the same norm. $K$, $DK$ are explicit
($DK_c = \partial_cP\,MP + P\,\partial_cM\,P + PM\,\partial_cP$), no AD through $p$.

*Factorised Jacobian.* $r$ depends on $W$ only through $(K, DK)$, which are
linear in $W$:
$\partial K_{ij}/\partial W_{(ar)f} = \phi_f P_{ia}P_{rj}$ and
$\partial DK_{ijc}/\partial W_{(ar)f} = \phi_f(\partial_cP_{ia}P_{rj} + P_{ia}\partial_cP_{rj}) + \partial_c\phi_f P_{ia}P_{rj}$.
So $\partial r/\partial W = G\cdot\partial(K,DK)/\partial W$ with
$G = \partial r/\partial(K,DK)$ from $N^2+N^3$ forward-mode tangents per point,
independent of the number of features: 5–16× cheaper than forward AD over $W$.

*Optimiser.* Adam warm-up on the first class only (2000 steps), then rounds of
dense Levenberg–Marquardt (Marquardt scaling, $\lambda\div 3$ on success,
$\times 4$ on failure, 8 iterations per round) on freshly resampled collocation
sets ($\ge 128$ points, $\ge 2\times$ as many residual rows as unknowns), with a
held-out estimate after every round. A class is accepted when the held-out
loss is $<10^{-20}$, or when LM stalls below $10^{-16}$; otherwise the ladder
climbs, warm-starting by copying coefficients monomial by monomial (never
re-running Adam on a warm start: its noise floor, $\sim 10^{-4}$, destroys an
LM optimum). Near these solutions LM converges only linearly (factor 2–10 per
round), consistent with degenerate zero sets (families of complex structures
and the diffeomorphism orbit), so the round budget matters more than the
iteration count.

### 5.4 Why a small *mean* Nijenhuis residual is not evidence on $S^6$, $S^4\times S^2$, $S^4\times S^1\times S^1$

$S^6\setminus\{pt\}\cong\mathbb C^3$, $S^2\times(S^4\setminus\{pt\})\cong\mathbb{CP}^1\times\mathbb C^2$
and $(S^4\setminus\{pt\})\times T^2\cong\mathbb C^2\times T^2$ are complex, and the
pulled-back structure extends as an ACS across a tube of radius $\varepsilon$
around the removed set (the obstructions lie in $\pi_3(\mathbb{CP}^3) =
\pi_5(\mathbb{CP}^3) = 0$). Inside the tube $|N|\sim 1/\varepsilon$ on a set of
measure $\sim\varepsilon^{c}$ ($c = 6$, resp. 4), so
$\mathbb E|N|^2\sim\varepsilon^{c-2}\to 0$ while $\sup|N|\to\infty$ —
whether or not a complex structure exists. A mean-square PINN is therefore
attracted to concentrated defects on exactly these cases (observed: the
conjugation run on $S^4\times S^2$ reached mean $|N| = 0.16$ with an adversarial
sup of $5\cdot10^{15}$, an affine conjugator gone singular at isolated points;
on $S^6$ 58 % of $\mathbb E|N|^2$ sits in 5 % of the points). Only a small and
resolution-stable **adversarial sup** of $|N|$ would be numerical evidence of a
CS there; §5.2's remark that a small residual is evidence of existence holds
for the sup, not for the mean.

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
  for class in [c, lin, bilin, quad_h, lin+c, lin+bilin, all2, (all2+trilin if >= 3 factors)]:
      for attempt in [plain, + annealed Dirichlet, + orientation penalty]   (next only if rejected):
          K = P poly_class(p) P ; minimise E||K^2 + P||^2  (3000 Adam steps)
          if adversarial_sup(||K^2+P||) < 1: accept; break
  K = P (poly + mlp)(p) P ; minimise E||K^2 + P||^2  (10k Adam + L-BFGS)
  if adversarial_sup(||K^2+P||) < 1:
      J* = K (-K^2 + I - P)^{-1/2}   (Newton–Schulz)      -> certified exact ACS
  else: NOT certified (S^4)
  # ---- CS unit ----------------------------------------------------------
  if no certified ACS: refuse
  cs-mode poly (default for new runs):
      for class in [deg 2, deg 3] (or [2t] with >= 3 factors):
          K = P M_W(p) P,  W = lsq(J*) on the first class, else copied from the previous class
          minimise E||N_K||^2 + 10 E||K^2+P||^2 : Adam (first class) + LM rounds (fresh points)
          accept if held-out loss < 1e-20 (or LM stalls below 1e-16)
      certify K, J = polar(K); report N of J (dense + adversarial)
  cs-mode conj (original):
      J = g2 g1 J* g1^-1 g2^-1,  g_s = exp(B_s) (or Cayley / I + B_s), B_s = P mlp_s(p) P
      minimise E ||N_J||^2  (Adam [+ tail sampling] + L-BFGS)
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
| orientation index $|o|$ of any exact ACS (§4.6) | $\ge 1$ (observed $1.000$–$1.012$ on every certified field) |
| degree of the Hopf ($S^3\times S^1$, $S^5\times S^1$) and Calabi–Eckmann ($S^3\times S^3$) structures as ambient fields $PMP$ | 2 (fit residual $10^{-15}$) |
| degree-2 floor of $\mathbb E\|N\|^2 + 10\,\mathbb E\|K^2+P\|^2$ on $S^6$, $S^4\times S^2$, $S^4\times S^1\times S^1$ (poly mode) | $\approx 4$–$11$ held out; the free field loses admissibility |

---

## 8b. Local campaign, 2026-09-30 (CPU, float64)

ACS units: seeds 7 and 8 on all ten cases. CS units: poly mode
(§5.3) from the unit's own certified ACS; on the three cases without a
near-integrable polynomial field, the conjugation mode (§5.1) with the
exponential map and tail-aware sampling as the sup-oriented best effort.
$N$ is always evaluated on an **exact** ACS (by construction, or the
polar-corrected certified field); $10^{-10}$ is the floor of the norm
evaluation ($\sqrt{\cdot + 10^{-20}}$).

| case | ACS: adversarial $\hat r$ (verdict) | CS: mode | mean $\|N\|$ | adversarial sup $\|N\|$ | chart CS: paper → now |
|---|---|---|---:|---:|---:|
| $S^2$ | $6\cdot10^{-6}$, $1\cdot10^{-5}$ (✓✓) | trivial ($d=2$) | $10^{-10}$ | $10^{-10}$ | $5.3\cdot10^{-4}\to10^{-10}$ |
| $S^1\times S^1$ | $1\cdot10^{-5}$, $6\cdot10^{-5}$ (✓✓) | trivial ($d=2$) | $10^{-10}$ | $10^{-10}$ | $1.9\cdot10^{-4}\to10^{-10}$ |
| $S^2\times S^2$ | $2.5\cdot10^{-4}$, $4.0\cdot10^{-4}$ (✓✓) | poly deg 2 (both seeds) | $1.0\cdot10^{-10}$ | $1.0\cdot10^{-10}$ | $4.0\cdot10^{-2}\to1.0\cdot10^{-10}$ |
| $S^3\times S^1$ | $1.1\cdot10^{-3}$, $8.4\cdot10^{-4}$ (✓✓) | poly deg 2 (both seeds) | $2.8$–$3.4\cdot10^{-10}$ | $0.7$–$1.2\cdot10^{-9}$ | $0.20\to8.4\cdot10^{-10}$ |
| $S^3\times S^3$ | $6.4\cdot10^{-4}$, $5.5\cdot10^{-4}$ (✓✓) | poly deg 2 (both seeds) | $1.0\cdot10^{-10}$ | $1.1\cdot10^{-10}$ | $7.63\to1.1\cdot10^{-10}$ |
| $S^5\times S^1$ | $3.1\cdot10^{-4}$, $5.9\cdot10^{-4}$ (✓✓) | poly deg 2 (both seeds) | $1.3$–$1.4\cdot10^{-10}$ | $2.0\cdot10^{-10}$ | $6.42\to2.7\cdot10^{-10}$ |
| $S^4$ | $2.00$, $2.00$ (✗✗: no ACS) | refused (no certified base) | — | — | $8.50\to$ — |
| $S^6$ | $4.4\cdot10^{-4}$, $4.1\cdot10^{-4}$ (✓✓) | conj exp + tail, 2500 steps (mean-square conj, 4000 steps) | $1.04$ ($0.38$) | $3.2$ ($6.9$) | $9.87\to2.8$ ($1.0$) |
| $S^4\times S^2$ | $2.7\cdot10^{-4}$, $3.1\cdot10^{-4}$ (✓✓) | conj exp + tail, 1000 of 2500 steps (affine conj, 4000) | $2.7$ ($0.16$) | $15$ ($5\cdot10^{15}$, singular $g$) | $6.28\to$ — |
| $S^4\times S^1\times S^1$ | $0.12$ (orientation attempt), $4.0\cdot10^{-3}$ (plain) (✓✓, both in `all2+trilin`) | conj exp + tail, 1000 of 2500 steps; poly `2t` not admissible | $1.5$ | $55$ | — |

On the four cases with a classical complex structure the search reaches the
floating-point floor in minutes, from a base ACS with $\|N\| = 1$–$17$. The
found structures are *not* round-orthogonal ($|J^{\top}J - P|\approx 0.07$–$0.19$,
singular values of $J$ in $[0.93, 1.07]$; on $S^3\times S^1$ and $S^3\times S^3$ the
defect is the same at every point), i.e. other members of the Hopf /
Calabi–Eckmann families than the specific structures of the offline check.
On $S^6$, $S^4\times S^2$ and $S^4\times S^1\times S^1$ no run comes near
integrability. Degree-2/3 polynomial classes hold no near-integrable field
(the free field loses admissibility), and the conjugation runs show the
trade-off of §5.4: the mean-square run on $S^6$ reaches mean $0.38$ with an
adversarial sup of $6.9$, the tail-aware run halves the sup ($3.2$) at mean
$1.04$. The two $S^4\times S^2$, $S^4\times S^1\times S^1$ tail runs were stopped at
1000 of 2500 steps because the machine was shared with another workload; the
collaborator branch `s4t2` reports RMS $0.29$ / max $0.78$ on
$S^4\times S^1\times S^1$ from a much longer HPC campaign. The complex
structures known on $S^6$ (Alpöge–Claude) and $S^2\times S^4$
(Liang–Liu–Zhang; arXiv:2609.05504) are built by algebro-geometric surgery
(Mumford degeneration of a modular family of 2-tori; blow-up, $\mu_2$-quotient,
Atiyah flop), so nothing suggests they are close to low-degree fields in the
round ambient picture; for $S^4\times T^2$ no existence result is known to us.

The conjugation mode, for comparison, stopped at $1.8\cdot10^{-4}$ ($S^2\times S^2$),
$1.3\cdot10^{-3}$ ($S^3\times S^1$), $4.6\cdot10^{-2}$ ($S^3\times S^3$) mean $\|N\|$ with the
same local budget.

## 8c. Degree-ladder measurement on exact ACS (2026-10-01)

Two exact-ACS measurement modes were added (`--cs-mode polar`, `--cs-mode conjpoly`),
both with an orthonormalised polynomial basis (SVD of the monomials over samples,
which removes the $|x_A| = 1$ redundancy; on $S^6$ the degree-$\le k$ basis has the
dimension of the harmonics of degree $\le k$) and the factorised exact Jacobian.

**Polar mode** ($J^* = \mathrm{polar}(K)$, $K = PM_WP$). The minimum of
$\mathbb E|N(J^*)|^2$ over the admissible part of a class is attained on its
boundary: on $S^6$ at degree 2, every LM step that lowers the Nijenhuis norm
does so by driving $K$ towards degeneracy between collocation points, and once
Newton–Schulz leaves its convergence region it returns a non-ACS whose
"Nijenhuis tensor" is artificially small. With exactness rows
$\sqrt\mu(J^{*2} + P)$, an admissibility barrier $\sqrt\mu\,\mathrm{relu}(|K^2+P| - 0.5)$,
tail points and dense rollback, the constrained floor stays at the starting
value ($|N|^2 \approx 384$): degree-2 free fields hold no admissible
improvement on $S^6$. Same pattern on $S^4\times S^2$.

**Polynomial-generator conjugation** ($J = e^{B}J_b e^{-B}$, $B = PM_WP$, $J_b$ the
unit's own certified ACS, $e^B$ by scaling and squaring). Exact ACS by construction,
$e^B$ always invertible, LM on resampled collocation sets ($\ge 25$ rows per
unknown) plus the 64 worst-$|N|$ points per round.

| $S^6$, degree of $B$ | coefficients | $\mathbb E|N|^2$ | mean $|N|$ | adversarial sup $|N|$ |
|---:|---:|---:|---:|---:|
| own ACS ($B = 0$) | — | 385 | 19.6 | 19.7 |
| 1 | 392 | 39 | 5.6 | 21 |
| 2 | 1715 | 4.4 | 1.9 | 12.6 |
| 3 | 5488 | 2.4 | 1.4 | 13.4 |

The mean decays (×9, then ×1.8) but the adversarial sup stops decreasing after
degree 2: additional degree buys mean at the cost of concentration (§5.4), not
an approach to a complex structure. For comparison, the collaborators' MLP
conjugation (≈200k parameters, first order, ≈160k steps on HPC) reports mean
0.067 and adversarial sup 1.5; no run, polynomial or MLP, has a sup that
decreases towards 0, which is the only signature that would count as evidence.

## 8d. Numerical reconstruction of the complex structures on $S^6$ and $S^2\times S^4$

**Status of the unbiased search.** No unsupervised, unbiased run on the round
ambient model reaches a complex structure on $S^6$ or $S^4\times S^2$ (§8b, §8c):
every ansatz has an adversarial Nijenhuis sup of order 1–10. What follows is a
*reconstruction*: the manifold is presented as the fibration of
[AC26]/[Engel] ($S^6$) and [GFL26] ($S^2\times S^4$), the function class and the
defining equations are taken from those papers, and a physics-informed solver
recovers the period map from residuals alone (no values of the period map or
of $J$ are given). The paper should present it as such.

**Model.** $X\supset (\mathfrak h_Z^\circ\times\mathbb R^4/\Lambda)/\Delta$ with
$\Delta = \Delta(3,4,\infty)$ acting by $(Z,x)\mapsto(g_jZ, A_jx)$,
$A_j = (T_j^{-1})^{t}$ [AC26, Rem. 3.20], and $J$ the x-invariant structure with
$(1,0)$-forms $dZ$ and $\Pi(Z)\,dx$, $\Pi = \begin{pmatrix}6\mu&\tau&1&0\\ \beta&\mu&0&1\end{pmatrix}$.
$S^2\times S^4$: the same period map on the index-two superlattice
$\Lambda^+ = \Lambda + \mathbb Z\hat\delta/2$ ($x = \mathrm{diag}(1,1,1,\tfrac12)x^+$, integral
monodromy $A_j^+$, equal to [GFL26, (3.2)]). Multiple fibres: log transforms on discs around
the fixed points $Z_1, Z_2$ (twisted translations leave the x-invariant $J$ unchanged).
Cusp: the toric quotient $Y_\epsilon/\Lambda$ of [AC26, §4] built from the same functions.

**Period-map PINN** (`extra/period_map.py`). In the coordinate $Z$ with the
distinguished cusp at $\infty$ and $g_0(Z) = Z-1$, the functions $\tau - Z$, $\mu$,
$\beta + Z$ are $g_0$-invariant and holomorphic on $\mathfrak h_Z$, hence power
series in $q = e^{2\pi i Z}$ on $|q|<1$; on the Ford domain $|q|\le 0.156$, so 21
coefficients suffice. Being exact q-series, they make $J$ integrable
identically; the loss is the residual of the six functional equations
(§3 of [AC26]) on discs around the fixed points **and on the gluing pairs of the
Ford domain** (collocating on the discs alone fits them and lets the high
coefficients blow up along the edges). $\tau$: Levenberg–Marquardt from
$\tau = Z + \text{const}$; $\mu$: linear least squares given $\tau$; $\beta$: linear
least squares given $\tau,\mu$, up to its free constant $c_0$ (the family
parameter; $c_0 = -3i$, so that $D<0$). Three seconds on a laptop.

**Certificate** (`extra/reconstruct_s6.py`, fresh points; coordinate norms in $(\mathrm{Re}Z, \mathrm{Im}Z, x)$):

| quantity | $S^6$ | $S^2\times S^4$ |
|---|---:|---:|
| $\tau(Z_1)-e^{i\pi/3}$, $\tau(Z_2)-i$, $\mu(Z_1)-\frac1{1+\rho}$, $\mu(Z_2)-\frac{1-i}2$ (not imposed) | $\le 5\cdot10^{-15}$ | same |
| functional-equation residual on the gluing pairs | $1.2\cdot10^{-13}$ | $9.9\cdot10^{-14}$ |
| $\|J^2+I\|$ | $2.5\cdot10^{-15}$ | $4.6\cdot10^{-15}$ |
| $\|N\|$ on the domain, adversarial sup | $1.1\cdot10^{-14}$ | $1.6\cdot10^{-14}$ |
| $\|N\|$ on the multiple-fibre discs | $\le 1.7\cdot10^{-14}$ | $\le 2.6\cdot10^{-14}$ |
| deck consistency $\|J(gZ) - D\Phi\,J(Z)\,D\Phi^{-1}\|$ | $\le 3.3\cdot10^{-12}$ | $\le 4.3\cdot10^{-12}$ |
| $\|N\|$ of the blended global field across the edges, adversarial | $\le 2.3\cdot10^{-12}$ | $\le 2.9\cdot10^{-12}$ |
| $\max D$ on the domain (every fibre a complex torus) | $-4.14$ | $-4.14$ |
| cusp: $C(q)$ single-valued; [AC26, Thm 4.5] $\epsilon_0 = \min(r/2, e^{-2\|R\|_\infty})$ | $1.8\cdot10^{-15}$; $1.8\cdot10^{-21} > 0$ | same |

so $\sup|N| \le 4.3\cdot10^{-12}$, six orders below the $10^{-6}$ target. Deep in
the cusp collar ($|q|$ down to $10^{-33}$) $|N|\le 2\cdot10^{-15}$. The other three
transport conventions for the deck action fail at $O(1)$–$O(100)$, so the lattice
action is pinned down numerically.

**What is and is not shown.** Shown: a complex structure, certified to
$4\cdot10^{-12}$, on the fibration model of $X$ (generic part, multiple-fibre
discs, gluing across all identifications), with the hypotheses of the toric
cusp filling verified for the learned period map (where the structure is the
standard toric one). Taken from the papers: that $X$ is diffeomorphic to $S^6$
(resp. $S^2\times S^4$), the topological twisting data of the fillings, and the
toric construction itself. Not shown: an explicit diffeomorphism to the round
$S^6$ — the reconstructed $J$ lives on the fibration model, which is why it is
invisible to the round-ambient search.

## 9. Files

| file | content |
|---|---|
| `src/geometry/ambient.py` | `AmbientManifold`: sampling, projector $P$ |
| `src/loss/ambient_loss.py` | ambient brackets, $N_{ij}$, $\|N\|$, $DJ$/$DP$ by forward AD |
| `src/geometry/explicit_acs.py` | case registry; classical ACS (offline checks / optional bases) |
| `extra/campaign_unit.py` | the unit: ladder, Dirichlet retry, MLP refinement, adversarial sup, Newton–Schulz polar, CS stage, chart protocol |
| `extra/campaign_summary.py` | before/after table with verdicts (poly-mode CS rows report the polar-corrected field) |
| `src/models/poly_field.py` | polynomial field $K = PMP$, explicit $DK$, intrinsic residual, factorised Jacobian |
| `src/optim/lm.py` | dense Levenberg–Marquardt |
| `extra/poly_cs.py` | the poly-mode CS stage (§5.3) |
| `extra/test_poly_cs.py` | tests: derivatives, Jacobian, LM, orientation lemma, three-factor class, conjugators |
| `extra/period_map.py` | period-map PINN for $\tau,\mu,\beta$ (q-series, functional-equation residuals, LM) |
| `extra/s6_reconstruction.py`, `extra/reconstruct_s6.py` | $J$ on the fibration model of $S^6$, $S^2\times S^4$ and its certificate |
| `extra/test_reconstruction.py` | tests of the reconstruction |
| `extra/fibration_cs.py`, `extra/atlas_cs.py` | exploratory: topology-only period search (plateau $2\cdot10^{-3}$), integrable atlas (not converged) |
| `scratch/submit_S4xS1xS1_{ACS,CS}.sh` | separate Lovelace submissions for the three-factor case |
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

[AC26] Claude and L. Alpöge, *The (3,4,∞) modular family of 2-tori, completed at its three special points, is a complex structure on S⁶*, https://alpo.ge/s6.pdf (2026). [Engel] P. Engel, *Complex structures on S⁶*, https://philip-engel.github.io/S6.pdf (2026). [GFL26] S. Guo, E. X. Fang, J. Lu, *A complex structure on S²×S⁴*, arXiv:2609.05504 (2026).
