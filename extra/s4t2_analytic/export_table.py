"""Write the classification of sphere products as a CSV table.

The classification itself is manuscript_draft/code/sphere_products.py (manuscript, Appendix A.1);
this script only imports it and tabulates its result. One row per product of spheres
S^{n_1} x ... x S^{n_k} of even dimension 2..12:

    dim, product, almost_complex, status, minimal_inputs, witness

status is one of
    no_acs       no almost complex structure (Sutherland / Borel-Serre)
    classical    complex structure known without any 2026 input
    needs_input  complex structure known, but every decomposition uses one of the
                 inputs [S6], [S2xS4], [LLZ], [S4xT2]
    open         almost complex, no known complex structure
minimal_inputs lists the inclusion-minimal input sets, separated by ';' (empty for
classical products); witness gives one decomposition into complex blocks for each.

Run:  python3 extra/s4t2_analytic/export_table.py [output.csv]
      (default output: extra/s4t2_analytic/data/sphere_products.csv)
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True      # leave no __pycache__ in manuscript_draft/
sys.path.insert(0, str(HERE.parent.parent / "manuscript_draft" / "code"))
from sphere_products import classify, name  # noqa: E402


def status(acs, decs):
    if not acs:
        return "no_acs"
    if not decs:
        return "open"
    return "classical" if frozenset() in decs else "needs_input"


def rows(max_dim=12):
    for dim in range(2, max_dim + 1, 2):
        for p, acs, decs in classify(dim):
            ordered = sorted(decs.items(), key=lambda kv: (len(kv[0]), sorted(kv[0])))
            yield {
                "dim": dim,
                "product": name(p),
                "almost_complex": acs,
                "status": status(acs, decs),
                "minimal_inputs": "; ".join("+".join(sorted(d)) for d, _ in ordered if d),
                "witness": "; ".join(" . ".join(w) for _, w in ordered),
            }


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "data" / "sphere_products.csv"
    table = list(rows())
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    print(f"wrote {len(table)} rows to {out}")


if __name__ == "__main__":
    main()
