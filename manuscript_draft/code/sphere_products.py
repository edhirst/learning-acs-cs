"""Status of complex structures on products of spheres S^{n_1} x ... x S^{n_k}.

Enumerates every product of spheres (standard smooth structures, n_i >= 1) of a
given even dimension and decides:

* whether it admits an almost complex structure (ACS):
    - all factors even: Sutherland 1965 (Thm 3.1, as cited in LLZ26): ACS iff the
      factors split into blocks S^2, S^6, S^2 x S^4.  (One factor: Borel-Serre.)
    - some factor odd: the product is parallelisable (Staples 1967), hence ACS.
* whether a complex structure is known, by searching for a decomposition into
  blocks that are known to be complex.  Each block records what it depends on.

Blocks (dependency tag in brackets):
    S^2                          CP^1                                     [classical]
    S^{2p+1} x S^{2q+1}          Calabi-Eckmann / Hopf / T^2              [classical]
    S^6                          Alpoge-Claude, Engel 2026                [S6]
    S^2 x S^4                    Guo-Fang-Lu 2026 or Liang-Liu-Zhang 2026  [S2xS4]
    S^4 x S^3 x S^{2q+1}         elliptic bundle over X_+ x P^q           [LLZ]
    S^4 x S^3 x S^4 x S^3        elliptic bundle over X_+ x X_+           [LLZ]
    S^4 x S^1 x S^1              our construction                         [S4xT2]

[LLZ] means: the complex S^2 x S^4 of Liang-Liu-Zhang together with its line
bundle L_+ (c_1 a generator).  It provides [S2xS4], and LLZ assume the S^6
threefold of Alpoge-Claude, so [LLZ] implies [S6] and [S2xS4].

Run:  python3 sphere_products.py        (prints dimensions 2..12, checks dim 8)
"""

from functools import lru_cache
from itertools import combinations

CLASSICAL, S6, S2S4, LLZ, S4T2 = "classical", "S6", "S2xS4", "LLZ", "S4xT2"
# LLZ's construction starts from the S^6 threefold and yields a complex S^2 x S^4,
# so [LLZ] implies [S6] and [S2xS4].
IMPLIES = {LLZ: {S6, S2S4}}


def partitions(n, largest=None):
    """Partitions of n into parts >= 1, as non-increasing tuples."""
    largest = n if largest is None else largest
    if n == 0:
        yield ()
        return
    for first in range(min(n, largest), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def remove(multiset, parts):
    """Remove the parts (a tuple) from the sorted tuple multiset, or return None."""
    m = list(multiset)
    for p in parts:
        if p not in m:
            return None
        m.remove(p)
    return tuple(m)


def candidate_blocks(ms):
    """Known complex blocks that fit inside the multiset ms, with their tag."""
    out = []
    if 2 in ms:
        out.append(((2,), CLASSICAL, "CP^1"))
    if 6 in ms:
        out.append(((6,), S6, "S^6"))
    if 2 in ms and 4 in ms:
        out.append(((4, 2), S2S4, "S^2xS^4"))
    odds = sorted({p for p in ms if p % 2})
    for a, b in combinations(odds + odds, 2):
        if a <= b and remove(ms, (a, b)) is not None:
            out.append(((a, b), CLASSICAL, f"CE(S^{a}xS^{b})"))
    if 4 in ms and 3 in ms:
        for c in odds:
            if remove(ms, (4, 3, c)) is not None:
                out.append(((4, 3, c), LLZ, f"S^4xS^3xS^{c}"))
        if remove(ms, (4, 4, 3, 3)) is not None:
            out.append(((4, 4, 3, 3), LLZ, "(S^4xS^3)^2"))
    if 4 in ms and remove(ms, (4, 1, 1)) is not None:
        out.append(((4, 1, 1), S4T2, "S^4xT^2"))
    # de-duplicate
    seen, uniq = set(), []
    for blk in out:
        key = (tuple(sorted(blk[0])), blk[1])
        if key not in seen:
            seen.add(key)
            uniq.append(blk)
    return uniq


def closure(deps):
    deps = set(deps)
    for d in list(deps):
        deps |= IMPLIES.get(d, set())
    deps.discard(CLASSICAL)
    return frozenset(deps)


def minimal(sets):
    """Inclusion-minimal members of a collection of frozensets."""
    sets = set(sets)
    return {s for s in sets if not any(t < s for t in sets)}


@lru_cache(maxsize=None)
def decompositions(ms):
    """Map: minimal dependency set -> one witness decomposition (list of block names)."""
    ms = tuple(sorted(ms))
    if not ms:
        return {frozenset(): ()}
    # every decomposition must use the largest part in some block
    top = ms[-1]
    results = {}
    for parts, tag, name in candidate_blocks(ms):
        if top not in parts:
            continue
        rest = remove(ms, tuple(sorted(parts)))
        if rest is None:
            continue
        for deps, witness in decompositions(rest).items():
            d = closure(set(deps) | {tag})
            if d not in results:
                results[d] = (name,) + witness
    keep = minimal(results)
    return {d: results[d] for d in keep}


def even_only_acs(ms):
    """Sutherland: product of even spheres is ACS iff it splits into S^2, S^6, S^2xS^4."""
    ms = list(ms)
    n4, n2 = ms.count(4), ms.count(2)
    others = [p for p in ms if p not in (2, 4, 6)]
    return not others and n4 <= n2


def has_acs(ms):
    if any(p % 2 for p in ms):
        return True  # parallelisable (Staples 1967)
    return even_only_acs(ms)


def name(ms):
    return " x ".join(f"S^{p}" for p in ms)


def classify(dim):
    rows = []
    for p in partitions(dim):
        ms = tuple(sorted(p))
        acs = has_acs(ms)
        decs = decompositions(ms) if acs else {}
        rows.append((p, acs, decs))
    return rows


def describe(decs):
    parts = []
    for d, w in sorted(decs.items(), key=lambda kv: (len(kv[0]), sorted(kv[0]))):
        tag = "classical" if not d else "+".join(sorted(d))
        parts.append(f"[{tag}] " + " . ".join(w))
    return "; ".join(parts)


def main(max_dim=12):
    for dim in range(2, max_dim + 1, 2):
        rows = classify(dim)
        no_acs = [r for r in rows if not r[1]]
        known = [r for r in rows if r[1] and r[2]]
        open_ = [r for r in rows if r[1] and not r[2]]
        print(f"\n=== dimension {dim}: {len(rows)} products ===")
        print(f"  no ACS ({len(no_acs)}): " + ", ".join(name(r[0]) for r in no_acs))
        needs_new = [r for r in known if frozenset() not in r[2]]
        if dim <= 8:
            print(f"  complex structure known ({len(known)}):")
            shown = known
        else:
            print(f"  complex structure known ({len(known)}); "
                  f"the {len(needs_new)} needing a 2026 input:")
            shown = needs_new
        for p, _, decs in shown:
            print(f"    {name(p):32s} {describe(decs)}")
        print(f"  OPEN: ACS exists, no known complex structure ({len(open_)}): "
              + (", ".join(name(r[0]) for r in open_) or "none"))

    # Assertions for dimension 8.
    rows = classify(8)
    assert len(rows) == 22
    no_acs = {r[0] for r in rows if not r[1]}
    assert no_acs == {(8,), (4, 4)}, no_acs
    assert all(r[2] for r in rows if r[1]), "an ACS product of dim 8 is still open"
    deps = {r[0]: r[2] for r in rows if r[1]}
    assert set(deps[(4, 3, 1)]) == {frozenset({LLZ, S6, S2S4})}
    assert set(deps[(4, 1, 1, 1, 1)]) == {frozenset({S4T2})}
    pre2026 = [r[0] for r in rows if r[1] and frozenset() not in r[2]]
    print("\nDimension 8 check passed: every product except S^8 and S^4xS^4 is complex,")
    print("given the inputs [S6], [S2xS4], [LLZ], [S4xT2].")
    print("  products needing a 2026 input: " + ", ".join(name(p) for p in pre2026))
    print("  S^4 x S^3 x S^1: only minimal input set is [LLZ] (which implies [S6], [S2xS4]);")
    print("  S^4 x T^4: only minimal input set is [S4xT2].")


if __name__ == "__main__":
    main()
