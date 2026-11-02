# learning-acs-cs

Code and run records for the paper **Learning Almost Complex and Complex Structures** by Edward Hirst, Tancredi Schettini Gherardini, Tomás S. R. Silva and Jakob Stein (2026).

The code searches numerically for almost complex and complex structures on spheres and products of spheres up to real dimension six, with neural networks and polynomial fields. The searches are not given any known solution; a separate supervised benchmark fits the known Hopf structures for comparison. The code also reconstructs complex structures on S⁶, S²×S⁴ and S⁴×T² from period data, and checks them numerically.

An almost complex structure is a field $J$ of endomorphisms of the tangent bundle with $J^2=-I$. It is a complex structure when its Nijenhuis tensor vanishes:

$$N_J(X,Y) = [JX,JY] - J[JX,Y] - J[X,JY] - [X,Y] = 0.$$

## Method

- **One global representation.** Each manifold is a product of unit spheres $M\subset\mathbb{R}^N$. A field is $K = P\,M(x)\,P$, where $P(x)$ projects onto the tangent space and $M(x)$ is an $N\times N$ matrix field: a polynomial in the ambient coordinates, a multilayer perceptron, or both. There are no charts, and $N_J$ is computed from ambient derivatives.
- **Almost complex stage.** Minimise $\mathbb{E}\,\lVert K^2+P\rVert^2$, trying polynomial feature classes of increasing size, then refine with a network. If the largest residual found by multi-start ascent is below 1, the field is corrected to an exact almost complex structure $J = K S^{-1/2}$ with $S = I-(K^2+P)$, computed by Newton–Schulz iteration.
- **Complex stage.** Starting from that structure $J_0$, minimise the Nijenhuis norm in one of two ways. A free polynomial field minimises $\mathbb{E}\,\lVert N_K\rVert^2 + 10\,\mathbb{E}\,\lVert K^2+P\rVert^2$ by Levenberg–Marquardt and is then corrected as above. Alternatively, $J = gJ_0g^{-1}$ with $g$ built from networks or from a polynomial generator, which keeps $J$ exactly almost complex.
- **Validation.** Means use fresh samples, and maxima come from multi-start gradient ascent of the pointwise residuals. These maxima are lower bounds on the true suprema, not proofs. Norms are evaluated as $(\lVert T\rVert^2+10^{-20})^{1/2}$, so $10^{-10}$ is a reporting floor.

## Results at a glance

| Manifold | Almost complex structure | Complex structure, direct search |
|---|---|---|
| S², T² | found by both seeds | found; automatic in dimension 2 |
| S⁴ | not found, as expected: S⁴ has none | not run, for lack of a starting structure |
| S²×S², S³×S¹, S³×S³, S⁵×S¹ | found by both seeds | found by both seeds, largest Nijenhuis residual at most 6.5×10⁻¹⁰ |
| S⁶, S⁴×S², S⁴×T² | found by both seeds | not found: the best candidate on each still has a largest Nijenhuis residual between 0.48 and 10.5 |

On S⁶, S²×S⁴ and S⁴×T², complex structures given by period data are reconstructed on fibration models (Section 5 of the paper). All recorded residuals are below 10⁻¹¹ (Table 6). They are measured in the coordinates of those models, so they are not comparable with the residuals of the direct search. These are reconstructions of given structures, not results of the search.

## Installation

```bash
git clone https://github.com/edhirst/learning-acs-cs.git
cd learning-acs-cs
python3 -m venv .venv
.venv/bin/pip install -r environment/requirements.txt pytest
source .venv/bin/activate
```

The environment must be in `.venv` at the repository root, because `experiments/run_all.py` starts every job with `.venv/bin/python`. Everything runs on CPU in double precision. The paper's records were made with Python 3.12.14, PyTorch 2.14.1, NumPy 2.5.3, SciPy 1.18.1 and Matplotlib 3.11.2 on a 10-core Apple M4 with 16 GB of memory.

The commands below are run from the repository root with the environment active.

## Tests

```bash
python -m pytest -q extra
```

This runs the 63 tests in about a minute.

## Rebuilding the paper's tables, figures and numbers

`outputs/paper_rerun` holds the records of every run in the paper ([description](outputs/paper_rerun/README.md)). These two commands rebuild the tables, figures and quoted numbers from them in about ten seconds, without training or model checkpoints:

```bash
python experiments/make_figures.py
python experiments/paper_numbers.py --json outputs/paper_rerun/numbers.json > outputs/paper_rerun/numbers.txt
```

The rebuilt tables and number files are identical to the committed ones, and the figure PDFs differ only in their embedded creation date. `numbers.txt` lists every number quoted in the paper with the record it comes from, keyed by table, figure or section:

```
T2.s2.s7.r_adv                                   4.198e-05   [acs_s2_s7: adv_acs_sup (after refinement)]
```

The paper's tables are condensed from the generated ones, which have more columns:

| Paper | Generated file in `outputs/paper_rerun` | Records in `outputs/paper_rerun/runs` |
|---|---|---|
| Figure 1, Table 2 | `figures/fig_acs_ladder.pdf`, `tables/tab_acs.tex` | `acs_<case>_s7`, `acs_<case>_s8` |
| Table 3 | `tables/tab_cs_poly.tex` | `cspoly_<case>_s<seed>` |
| Figure 2 | `figures/fig_cs_conj_vs_poly.pdf` | `csconjexp_*`, `csconjaff_*`, `cspoly_*_s7` |
| Table 4 | `tables/tab_supervised.tex` | `hopf_*`, `hopf_class`, `csconj*_s3xs1`, `csconj*_s5xs1` |
| Table 5 | `tables/tab_hard.tex` | the seed-7 runs on `s6`, `s4xs2` and `s4xs1xs1` |
| Figure 3 | `figures/fig_defect.pdf` | `cs_analysis` |
| Figure 4 | `figures/fig_recon_convergence.pdf` | `reconstruction_K` |
| Table 6 | `tables/tab_recon.tex`, `tables/tab_s4t2_recon.tex` | `reconstruction`, `s4t2_reconstruction` |
| Figure 5 | `figures/fig_s4t2_recon.pdf` | `s4t2_reconstruction` |

`make_figures.py` also draws some figures and tables that the paper does not use.

## Rerunning the experiments

```bash
PAPER_OUT=outputs/rerun python experiments/run_all.py --budget 9
PAPER_OUT=outputs/rerun python experiments/make_figures.py
```

`run_all.py` runs the 77 jobs of the campaign in dependency order, keeping at most `--budget` threads busy, and writes to `$PAPER_OUT` (default `outputs/paper_rerun`, the paper's records). A job counts as done once its output exists, so repeating the command resumes an interrupted campaign. `--dry` lists the jobs, and `--only` selects them by name prefix, for example `--only acs_s2xs2,cspoly_s2xs2`.

The paper's campaign took about 26 hours with `--budget 9` on the machine above; the four supervised fits on S⁵×S¹ take about 19 hours each. Runs are seeded and `run_all.py` fixes the thread count of each job, so on the same software and hardware a job reproduces its record exactly (checked for `acs_s2xs2_s7` and `cspoly_s2xs2_s7`). Another thread count, library version or machine can change the numbers slightly.

| Jobs (run folders) | Experiment |
|---|---|
| `acs_<case>_s<seed>` | Almost complex stage; seeds 7 and 8, plus 9 and 10 on S³×S¹ |
| `cspoly_<case>_s<seed>` | Complex stage with a free polynomial field and Levenberg–Marquardt |
| `csconjexp_*`, `csconjaff_*`, `csconjtail_*` | Complex stage by network conjugation: exponential, affine, and exponential with tail-aware sampling |
| `conjpoly_*` | Complex stage by conjugation with a polynomial generator |
| `hopf_*` | Supervised fits of the Hopf structures on S³×S¹ and S⁵×S¹ |
| `reconstruction`, `reconstruction_K` | Reconstruction on S⁶ and S²×S⁴, and the resolution study of the period solve |
| `s4t2_reconstruction` | Reconstruction of the S⁴×T² period model |
| `s4t2inv*` | Translation-invariant searches on S⁴×T² |
| `offline_checks`, `hopf_class`, `cs_analysis`, `figures` | Checks, analysis of the trained fields, then figures and tables |

## Running one search

```bash
python extra/campaign_unit.py --case s2xs2 --task acs --seed 7 --out outputs/my_runs/acs_s2xs2_s7
python extra/campaign_unit.py --case s2xs2 --task cs --cs-mode poly --poly-degrees 2,3 --seed 7 \
    --base-from outputs/my_runs/acs_s2xs2_s7 --out outputs/my_runs/cspoly_s2xs2_s7
```

The first command finds an almost complex structure on S²×S², and the second a complex structure starting from it; together they take about four minutes. They use every available thread, so their numbers differ slightly from the records, which `run_all.py` made with one thread per job. The cases are `s2`, `s4`, `s6`, `s1xs1` (T²), `s2xs2`, `s3xs1`, `s3xs3`, `s4xs2`, `s5xs1` and `s4xs1xs1` (S⁴×T²). Each run writes `result.json`, `progress.jsonl` and model checkpoints to its folder. `--help` lists the options, and `experiments/run_all.py` gives the settings of every job in the paper.

## Repository layout

| Path | Contents |
|---|---|
| `extra/campaign_unit.py` | One search: a (manifold, task, seed) unit with the almost complex and complex stages |
| `extra/poly_cs.py` | Polynomial complex-structure modes (free field and polynomial generator), fitted by Levenberg–Marquardt |
| `extra/ambient_supervised_hopf.py` | Supervised fits of the Hopf structures |
| `extra/hopf_invariant.py` | Homotopy class of almost complex structures on S³×S¹ relative to the Hopf structure |
| `extra/period_map.py`, `extra/s6_reconstruction.py`, `extra/reconstruct_s6.py` | Period-map solve and reconstruction on S⁶ and S²×S⁴ |
| `extra/s4t2_reconstruction.py`, `extra/s4t2_invariant.py` | S⁴×T²: reconstruction of the period model, and the translation-invariant reduction to S⁴ |
| `extra/test_*.py` | Tests |
| `src/` | Geometry of sphere products, Nijenhuis tensors, polynomial fields, Levenberg–Marquardt |
| `experiments/` | Scheduler, figures and tables, number index, analysis and S⁴×T² drivers |
| `outputs/paper_rerun/` | Records, tables and figures of the paper's campaign |
| `docs/` | Development notes on the mathematics of the ambient method and on the supervised Hopf fits |

`train.py`, `evaluate.py`, `configs/`, `extra/good_models/` and `analysis_v1.ipynb` belong to the earlier chart-based network mentioned in Remark 1.1 of the paper; `extra/test_all.py` still tests it. `pbs/` and `jobs_*.txt` are cluster scripts and job lists from earlier campaigns. The other scripts in `extra/` are tools from the development of the method. None of these produced the paper's results.

## Citation

```bibtex
@misc{HirstSchettiniGherardiniSilvaStein2026,
  author = {Hirst, Edward and Schettini Gherardini, Tancredi and Silva, Tom{\'a}s S. R. and Stein, Jakob},
  title  = {Learning Almost Complex and Complex Structures},
  year   = {2026}
}
```

## License

MIT; see [LICENSE](LICENSE).
