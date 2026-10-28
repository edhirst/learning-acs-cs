"""Homotopy class of an ACS on S³×S¹ relative to the Hopf structure J_t (z ↦ e^{2πc} z).

The exact-ACS conjugation J = e^{B₂} e^{B₁} J_base e^{−B₁} e^{−B₂} is homotopic to J_base, so a
CS search or a supervised fit by conjugation only reaches structures in the base point's class.
On S³×S¹ the first Chern class carries no information (H² = 0, cf. extra/chern_class_ambient.py);
the classes differ on the 3-cell, by π₃(S²) = ℤ:

In the global frame V_k = |x| e_k of (ℂ²∖0)/⟨z ↦ e^{2πc} z⟩, with ambient components
(e_k − q_k q, (q_k/c) ∂_θ), J_t is the constant matrix I. An ACS J of the same orientation becomes
M(p) = V⁺ J V, a complex structure on ℝ⁴; its orthogonalisation (with respect to the J-Hermitian
metric 1 + MᵀM) is a unit anti-self-dual 2-form, i.e. a point f(p) ∈ S², and J_t is a constant
map. The Hopf invariant of f on a slice S³×{θ} — the linking number of two preimage curves,
computed on a Kuhn triangulation of S³ — is a homotopy invariant of J relative to J_t: zero is
necessary for J to be conjugated to J_t.

Validation (python extra/hopf_invariant.py): the Hopf map (+1), its precomposition with the
quaternion conjugation (−1) and squaring (+2), a null-homotopic map (0); the pull-back of J_t by
the orientation-preserving isometry (q, θ) ↦ (q̄, −θ), an exact complex structure, has +1.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import os
import sys
import math
import itertools

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

C_MAIN = 1.0 / (2 * math.pi)


def frame(p, c=C_MAIN):
    """(B, 6, 4): columns V_k = (e_k − q_k q, (q_k/c) ∂_θ) at ambient points p = (q, w) of S³×S¹."""
    q, w = p[:, :4], p[:, 4:]
    tang = torch.eye(4, dtype=p.dtype).unsqueeze(0) - q.unsqueeze(2) * q.unsqueeze(1)
    R = torch.stack([-w[:, 1], w[:, 0]], dim=1)
    return torch.cat([tang, (R.unsqueeze(1) * (q / c).unsqueeze(2)).transpose(1, 2)], dim=1)


def s2_components(J, p, c=C_MAIN):
    """Self-dual and anti-self-dual parts (each (B, 3)) of the orthogonalised frame matrix of J."""
    V = frame(p, c)
    M = torch.linalg.pinv(V) @ J @ V
    g = torch.eye(4, dtype=M.dtype) + M.transpose(1, 2) @ M          # J-Hermitian metric
    ev, U = torch.linalg.eigh(g)
    Jo = (U * ev.sqrt().unsqueeze(1)) @ U.transpose(1, 2) @ M @ (U / ev.sqrt().unsqueeze(1)) @ U.transpose(1, 2)
    sd = torch.stack([Jo[:, 0, 1] + Jo[:, 2, 3], Jo[:, 0, 2] - Jo[:, 1, 3], Jo[:, 0, 3] + Jo[:, 1, 2]], 1) / 2
    asd = torch.stack([Jo[:, 0, 1] - Jo[:, 2, 3], Jo[:, 0, 2] + Jo[:, 1, 3], Jo[:, 0, 3] - Jo[:, 1, 2]], 1) / 2
    return sd, asd


def sphere_mesh(n):
    """Kuhn triangulation of the boundary of [−1, 1]⁴ (n cells per edge), radially projected to S³.
    The Kuhn subdivision restricts consistently to shared 2-faces: the tets form a closed 3-manifold."""
    tets = []
    for axis in range(4):
        free = [a for a in range(4) if a != axis]
        cells = np.array(list(itertools.product(range(n), repeat=3)))
        for side in (0, n):
            base = np.zeros((len(cells), 4), dtype=np.int64)
            base[:, axis] = side
            base[:, free] = cells
            for perm in itertools.permutations(free):
                chain = [base.copy()]
                for a in perm:
                    nxt = chain[-1].copy()
                    nxt[:, a] += 1
                    chain.append(nxt)
                tets.append(np.stack(chain, 1))                         # (cells, 4 vertices, 4 coords)
    T = np.concatenate(tets)
    lin = ((T[..., 0] * (n + 1) + T[..., 1]) * (n + 1) + T[..., 2]) * (n + 1) + T[..., 3]
    uniq, inv = np.unique(lin.ravel(), return_inverse=True)
    X = np.stack(np.unravel_index(uniq, (n + 1,) * 4), 1) * (2.0 / n) - 1.0
    return X / np.linalg.norm(X, axis=1, keepdims=True), inv.reshape(T.shape[:2])


def preimage_segments(X, f, tets, y):
    """Oriented polygonal preimage f⁻¹(y): affine interpolation, in each tet, of the stereographic
    chart of S² centred at y. Returns segment start and end points on S³ and the number of tets
    crossed an odd number of times (0 for a well-resolved map)."""
    e1 = np.cross(y, [1.0, 0, 0] if abs(y[0]) < 0.9 else [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(y, e1)
    den = 1 + f @ y
    v = np.stack([f @ e1, f @ e2], 1) / np.maximum(den, 1e-12)[:, None]
    T = tets[~(den[tets] < 0.2).any(1)]                                 # the chart is unreliable near −y
    vt, xt = v[T], X[T]
    hit, pts = [], []
    for i in range(4):                                                  # face opposite vertex i
        o = [j for j in range(4) if j != i]
        Mx = np.concatenate([vt[:, o, :].transpose(0, 2, 1), np.ones((len(T), 1, 3))], axis=1)
        ok = np.abs(np.linalg.det(Mx)) > 1e-14
        lam = np.zeros((len(T), 3))
        lam[ok] = np.linalg.solve(Mx[ok], np.broadcast_to([0.0, 0.0, 1.0], (int(ok.sum()), 3))[..., None])[..., 0]
        hit.append(ok & (lam >= 0).all(1))
        pts.append(np.einsum('tj,tjk->tk', lam, xt[:, o, :]))
    hit = np.stack(hit, 1)
    A, B = [], []
    for k in np.nonzero(hit.sum(1) == 2)[0]:
        i, j = np.nonzero(hit[k])[0]
        a, b = pts[i][k], pts[j][k]
        E = (xt[k, 1:] - xt[k, 0]).T                                     # (4, 3) edges of the tet
        D = (vt[k, 1:] - vt[k, 0]).T                                     # (2, 3) chart increments
        # Pontryagin orientation: (t, a, b) positive on S³ when (df a, df b) is positive on S²
        t = E @ np.cross(D[0], D[1]) * np.sign(np.linalg.det(np.column_stack([xt[k].mean(0), E])))
        forward = (b - a) @ t > 0
        A.append(a if forward else b)
        B.append(b if forward else a)
    return np.array(A), np.array(B), int((hit.sum(1) % 2 == 1).sum())


def linking(A1, B1, A2, B2):
    """Gauss linking number of two closed polygonal links in S³: exact segment-pair solid angles
    (Klenin & Langowski 2000) after an orientation-preserving stereographic projection from a
    pole away from both links."""
    allp = np.concatenate([A1, B1, A2, B2])
    cand = np.random.default_rng(0).normal(size=(4000, 4))
    cand /= np.linalg.norm(cand, axis=1, keepdims=True)
    pole = cand[np.argmin((allp @ cand.T).max(0))]
    Q, _ = np.linalg.qr(np.column_stack([pole, np.eye(4)[:, :3]]))
    Bs = Q[:, 1:].copy()
    if np.linalg.det(np.column_stack([pole, Bs])) > 0:    # pole e₄ with e₁..e₃: det < 0 preserves orientation
        Bs[:, 0] *= -1
    st = lambda x: (x @ Bs) / (1 - x @ pole)[:, None]
    a1, b1, a2, b2 = st(A1), st(B1), st(A2)[None], st(B2)[None]
    unit = lambda u: u / np.linalg.norm(u, axis=-1, keepdims=True)
    asin = lambda u, w: np.arcsin(np.clip((u * w).sum(-1), -1, 1))
    total = 0.0
    for s in range(0, len(a1), 400):
        p1, p2 = a1[s:s + 400, None], b1[s:s + 400, None]
        r13, r14, r23, r24 = a2 - p1, b2 - p1, a2 - p2, b2 - p2
        n1, n2 = unit(np.cross(r13, r14)), unit(np.cross(r14, r24))
        n3, n4 = unit(np.cross(r24, r23)), unit(np.cross(r23, r13))
        om = asin(n1, n2) + asin(n2, n3) + asin(n3, n4) + asin(n4, n1)
        total += np.nansum(om * np.sign((np.cross(b2 - a2, p2 - p1) * r13).sum(-1)))
    return total / (4 * np.pi)


def hopf_invariant(X, f, tets, trials=3, seed=1):
    """Linking numbers of the preimages of `trials` random pairs of values (all equal, and integers,
    for a well-resolved map); f: (V, 3) unit vectors at the mesh vertices X."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(trials):
        ys = rng.normal(size=(2, 3))
        ys /= np.linalg.norm(ys, axis=1, keepdims=True)
        A1, B1, odd1 = preimage_segments(X, f, tets, ys[0])
        A2, B2, odd2 = preimage_segments(X, f, tets, ys[1])
        out.append({'linking': float(linking(A1, B1, A2, B2)) if len(A1) and len(A2) else 0.0,
                    'segments': [len(A1), len(A2)], 'odd_tets': odd1 + odd2})
    return out


def relative_class(J_fn, Jt_fn, theta=0.3, n=24, mesh=None, chunk=4096):
    """Hopf invariant of J relative to J_t on the slice S³×{θ}, for J of the orientation of J_t
    (anything else raises). Returns (invariant or None if the trials disagree, details)."""
    X, tets = mesh if mesh is not None else sphere_mesh(n)
    P = torch.tensor(np.c_[X, np.tile([math.cos(theta), math.sin(theta)], (len(X), 1))], dtype=torch.float64)
    with torch.no_grad():
        sdt, asdt = s2_components(Jt_fn(P[:1]), P[:1])
        use_sd = bool(sdt[0].norm() > asdt[0].norm())
        ft = (sdt if use_sd else asdt)[0].numpy()
        sd, asd = map(torch.cat, zip(*[s2_components(J_fn(P[i:i + chunk]), P[i:i + chunk])
                                       for i in range(0, len(P), chunk)]))
    f = (sd if use_sd else asd).numpy()
    if not (np.linalg.norm(f, axis=1) > 0.5).all():
        raise ValueError('J does not have the orientation of the target')
    f /= np.linalg.norm(f, axis=1, keepdims=True)
    ang = np.degrees(np.arccos(np.clip(f @ (ft / np.linalg.norm(ft)), -1, 1)))
    trials = hopf_invariant(X, f, tets)
    lk = [t['linking'] for t in trials]
    inv = int(round(lk[0])) if max(lk) - min(lk) < 0.1 and abs(lk[0] - round(lk[0])) < 0.1 else None
    return inv, {'theta': theta, 'angle_to_target_deg': [float(ang.min()), float(ang.max())], 'trials': trials}


if __name__ == '__main__':
    from src.geometry.hopf_ambient import hopf_ambient_acs
    Xm, tm = sphere_mesh(24)

    def hopf(q):
        z1, z2 = q[:, 0] + 1j * q[:, 1], q[:, 2] + 1j * q[:, 3]
        w = z1 * np.conj(z2)
        return np.stack([np.abs(z1) ** 2 - np.abs(z2) ** 2, 2 * w.real, 2 * w.imag], 1)
    qsq = lambda q: np.stack([q[:, 0] ** 2 - (q[:, 1:] ** 2).sum(1), *(2 * q[:, 0:1] * q[:, 1:]).T], 1)
    for name, fm, h in (('Hopf map', hopf(Xm), 1), ('Hopf ∘ conjugation', hopf(Xm * [1, -1, -1, -1]), -1),
                        ('Hopf ∘ squaring', hopf(qsq(Xm)), 2), ('null-homotopic', np.c_[Xm[:, :2], 2 + Xm[:, 2]], 0)):
        fm = fm / np.linalg.norm(fm, axis=1, keepdims=True)
        lk = [round(t['linking'], 3) for t in hopf_invariant(Xm, fm, tm, trials=2)]
        print(f"{name:22s} expected {h:+d}: linking {lk}")
    Jt = lambda p: hopf_ambient_acs(p, 2, C_MAIN)
    S = torch.diag(torch.tensor([1., -1, -1, -1, 1, -1], dtype=torch.float64))
    for name, J, h in (('J_t', Jt, 0), ('−J_t', lambda p: -Jt(p), 0),
                       ('(q̄, −θ)*J_t', lambda p: S @ Jt(p @ S) @ S, 1)):
        inv, info = relative_class(J, Jt, mesh=(Xm, tm))
        print(f"{name:22s} expected {h:+d}: relative Hopf invariant {inv}")
