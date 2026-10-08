# Environment for the analytic-construction checks

The checks of the manuscript (`manuscript_draft/code/`) and the tools in `extra/s4t2_analytic/` need only
Python and `sympy`. They run on a laptop in about a second and need no GPU.

The symbolic campaign (`extra/s4t2_explicit/`) needs more, including numpy and torch. If you already have
that environment, it can run these checks too, as long as it has sympy 1.14.

## Create the virtual environment

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r extra/s4t2_analytic/environment/requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\activate` instead. `.venv/` is excluded by
`.gitignore`.

## Dependencies

| Package | Version | Used for |
|---|---|---|
| Python | ≥ 3.9 (tested with 3.10.13) | sympy 1.14 needs 3.9 or later |
| `sympy` | 1.14.0 | Exact symbolic algebra in `verify_construction.py` |
| `mpmath` | 1.3.0 | Required by sympy |

`sphere_products.py` and `export_table.py` use only the Python standard library.

The versions are pinned so that the outputs match `data/` exactly; the ledger baseline also records
sympy 1.14. Any recent sympy should still pass every check, but its printed expressions may be formatted
differently.

## Check the installation

```bash
python -c "import sympy; print(sympy.__version__)"
```

This should print `1.14.0`.
