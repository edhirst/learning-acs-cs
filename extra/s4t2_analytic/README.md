# `extra/s4t2_analytic/`: code and data for the analytic construction

This folder supports the closed-form construction of the complex threefold Y ≅ S⁴ × T² in
*An Explicit Complex Structure on S⁴ × T²*. It holds the environment, a runner and the reference outputs
for the manuscript's two verification scripts.

**The full write-up of the construction and its proofs is
[`docs/s4t2_analytic_construction.md`](../../docs/s4t2_analytic_construction.md).** Its Appendix A maps
every check to the statement it verifies. The manuscript itself is `manuscript_draft/main.tex`.

The scripts are **not copied** here. They live in `manuscript_draft/code/`, which is their single source of
truth (`AGENTS.md` §3.3), and everything in this folder calls them there.

## Layout

| Path | Contents |
|---|---|
| `manuscript_draft/code/verify_construction.py` | 15 exact `sympy` checks of the construction, matching manuscript Appendix A. They cover the lattice determinant, monodromy, the chart at ∞, the twisting section, the T²-action, the canonical bundle, the toric charts and dual bases, the local determinant identity, the exact heights P₀ and P_∞, the periods and Néron–Severi group, the Kodaira–Spencer derivative, and the admissible deformation region. All arithmetic is exact, with no floating-point tolerance. |
| `manuscript_draft/code/sphere_products.py` | The products of spheres of even dimension 2–12: it decides almost complexity (Sutherland, Staples), searches decompositions into known complex blocks, records the minimal inputs each product needs, and asserts the dimension-8 classification (Corollary D) |
| `export_table.py` | Imports `sphere_products.py` and writes its classification as a CSV table |
| `run_all.sh` | Runs both scripts and the export, and regenerates `data/` |
| `data/verify_construction.txt` | Reference output of the 15 checks |
| `data/sphere_products.txt` | Reference output of the classification, by dimension |
| `data/sphere_products.csv` | Classification table, one row per product (159 rows), with columns `dim, product, almost_complex, status, minimal_inputs, witness` |
| `environment/` | Virtual-environment setup (`README.md`, `requirements.txt`) |

## Running the code

Set up the environment as described in [`environment/README.md`](environment/README.md):

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r extra/s4t2_analytic/environment/requirements.txt
```

Then, from the repository root, run everything:

```bash
extra/s4t2_analytic/run_all.sh
```

Or run each script on its own:

```bash
python manuscript_draft/code/verify_construction.py       # ends with "all checks passed"
(cd manuscript_draft/code && python sphere_products.py)   # ends with "Dimension 8 check passed ..."
python extra/s4t2_analytic/export_table.py                # writes data/sphere_products.csv
```

These are the two checks of `AGENTS.md` §3.3. The LaTeX build is the third; it is run from
`manuscript_draft/` and is independent of this folder.

### How the code works

**`verify_construction.py`** runs from top to bottom.
- Each numbered block builds the relevant objects as `sympy` expressions in real coordinates
  t = a + ib: the lattice vectors, the transformed generators, the filling data.
- It then `assert`s the identity stated in the manuscript and prints one summary line.
- Every check is an exact identity, apart from a few inequalities, and those are also decided exactly:
  - check 1 reduces the determinant to 2π(1+r²)(2r² − (1−r²) ln r); its positivity for all r > 0 is a
    sign argument made in the manuscript, which the script spot-checks at r = 1, ½, 2;
  - checks 11 and 14 compare exact values with rational bounds at the stated points.
- The script ends with `all checks passed`. A failed assertion stops it with a traceback.

**`sphere_products.py`** works in these steps.
- It generates the partitions of each dimension.
- It decides almost complexity: Sutherland's criterion when all factors are even, parallelisability when
  some factor is odd.
- By memoised recursion it finds every way to split each product into complex blocks: ℂP¹, Calabi–Eckmann
  S^{2p+1} × S^{2q+1}, S⁶, S² × S⁴, S⁴ × S³ × S^{2q+1}, (S⁴ × S³)² and S⁴ × T².
- For each product it keeps the inclusion-minimal sets of inputs: the 2026 constructions [S6], [S2xS4],
  [LLZ], and the Main Theorem [S4xT2].
- Finally it asserts the dimension-8 table of Corollary D.

Each run takes about a second. The outputs should match the files in `data/` exactly. If a script in
`manuscript_draft/code/` changes, rerun `run_all.sh` and commit the updated `data/` with it.
