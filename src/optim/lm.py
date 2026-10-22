"""
Dense Levenberg–Marquardt for small parameter tensors (≲ 10^4 entries).

Minimises  L(θ) = Σ_c w_c |r_c(θ)|²  over a list of cached batches c. Exact
Jacobians J_c = ∂r_c/∂θ by forward-mode AD (vmap over chunks of the identity
basis), normal equations accumulated chunk by chunk, Marquardt scaling
D = diag(JᵀJ), step accepted iff L decreases (λ ÷ 3), otherwise λ × 4 and the
step is recomputed from the same Jacobian. For zero-residual problems (an
exact solution in the class) this converges quadratically, where first-order
and quasi-Newton methods stall on the ill-conditioned valley floor.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import time
from typing import Callable, List, Optional, Sequence

import torch


def jacobian_fwd(f: Callable[[torch.Tensor], torch.Tensor], theta: torch.Tensor,
                 tchunk: int = 128) -> torch.Tensor:
    """(m, n) Jacobian of f: R^n -> R^m at theta (flat) by chunked forward mode."""
    n = theta.numel()
    cols = []
    for s in range(0, n, tchunk):
        e = min(n, s + tchunk)
        basis = torch.zeros(e - s, n, dtype=theta.dtype, device=theta.device)
        basis[torch.arange(e - s), torch.arange(s, e)] = 1.0
        jv = torch.func.vmap(lambda v: torch.func.jvp(f, (theta,), (v,))[1])(basis)
        cols.append(jv)
    return torch.cat(cols, dim=0).T


def lm_minimise(residual_fn: Callable[[torch.Tensor, object], torch.Tensor],
                theta0: torch.Tensor, caches: Sequence[object],
                weights: Optional[Sequence[float]] = None, n_iter: int = 30,
                lam0: float = 1e-3, tchunk: int = 128, tol: float = 0.0,
                time_limit_s: Optional[float] = None,
                callback: Optional[Callable[[int, torch.Tensor, dict], None]] = None,
                jac_fn: Optional[Callable[[torch.Tensor, object], tuple]] = None):
    """Returns (theta, log). residual_fn(theta_shaped, cache) -> flat residual.
    `weights` default to 1/len(caches) (i.e. the mean over batches).
    `jac_fn(theta_shaped, cache) -> (r, J)` (optional): an exact Jacobian
    supplied by the model (e.g. the factorised polynomial-field Jacobian);
    default = chunked forward-mode AD of residual_fn."""
    shape = theta0.shape
    theta = theta0.detach().clone().reshape(-1)
    n = theta.numel()
    weights = list(weights) if weights is not None else [1.0 / len(caches)] * len(caches)

    def total(th):
        with torch.no_grad():
            return sum(wc * (residual_fn(th.view(shape), c) ** 2).sum().item()
                       for c, wc in zip(caches, weights))

    L = total(theta)
    lam = lam0
    log: List[dict] = [{'iter': 0, 'loss': L, 'lam': lam}]
    t0 = time.time()
    for it in range(1, n_iter + 1):
        JTJ = torch.zeros(n, n, dtype=theta.dtype)
        g = torch.zeros(n, dtype=theta.dtype)
        for c, wc in zip(caches, weights):
            if jac_fn is not None:
                r, J = jac_fn(theta.view(shape), c)
            else:
                f = lambda th, c=c: residual_fn(th.view(shape), c)
                with torch.no_grad():
                    r = f(theta)
                J = jacobian_fwd(f, theta, tchunk)
            JTJ += wc * (J.T @ J)
            g += wc * (J.T @ r)
            del J
        diag = torch.diagonal(JTJ).clone()
        D = diag.clamp_min(1e-10 * diag.mean().item() + 1e-300)
        accepted = False
        delta = torch.zeros_like(theta)          # no step if lam is already saturated
        while lam < 1e12:
            A = JTJ + lam * torch.diag(D)
            try:
                Lc = torch.linalg.cholesky(A)
                delta = -torch.cholesky_solve(g.unsqueeze(1), Lc).squeeze(1)
            except RuntimeError:
                delta = -torch.linalg.lstsq(A, g.unsqueeze(1)).solution.squeeze(1)
            L_new = total(theta + delta)
            if L_new < L:
                theta = theta + delta
                L, lam, accepted = L_new, max(lam / 3.0, 1e-15), True
                break
            lam *= 4.0
        rec = {'iter': it, 'loss': L, 'lam': lam, 'accepted': accepted,
               'step_norm': delta.norm().item(), 'grad_norm': g.norm().item(),
               'elapsed_s': time.time() - t0}
        log.append(rec)
        if callback is not None:
            callback(it, theta.view(shape), rec)
        if not accepted or L <= tol:
            break
        if time_limit_s is not None and time.time() - t0 > time_limit_s:
            break
    return theta.view(shape), log


def lm_cg_minimise(residual_fn, factor_fn, theta0, caches, weights=None, n_iter=10,
                   lam0=1e-3, cg_iters=200, cg_tol=1e-10, callback=None):
    """Matrix-free Levenberg-Marquardt: per iteration, build one factorised
    Jacobian per cache (`factor_fn(theta_shaped, cache)` -> object with
    residual(), matvec(v), rmatvec(u)), then solve
        (J^T J + lam * (mean diag proxy) I) delta = -J^T r
    by conjugate gradients (Levenberg damping, scaled by the mean of a
    randomised diag(J^T J) estimate). Accept if the true loss drops."""
    shape = theta0.shape
    theta = theta0.detach().clone().reshape(-1)
    weights = list(weights) if weights is not None else [1.0 / len(caches)] * len(caches)

    def total(th):
        with torch.no_grad():
            return sum(wc * (residual_fn(th.view(shape), c) ** 2).sum().item()
                       for c, wc in zip(caches, weights))
    L, lam, log = total(theta), lam0, []
    log.append({'iter': 0, 'loss': L, 'lam': lam})
    t0 = time.time()
    for it in range(1, n_iter + 1):
        facs = [factor_fn(theta.view(shape), c) for c in caches]
        sw = [w ** 0.5 for w in weights]
        g = sum(s * f.rmatvec(s * f.residual()) for f, s in zip(facs, sw))
        # scale of J^T J: Hutchinson estimate of its mean diagonal
        z = torch.randn(4, theta.numel(), dtype=theta.dtype)
        dscale = sum((sum((s * f.matvec(zz)) .pow(2).sum() for f, s in zip(facs, sw))).item()
                     for zz in z) / (4 * theta.numel())
        accepted = False
        delta = torch.zeros_like(theta)          # no step if lam is already saturated
        while lam < 1e12:
            mu = lam * dscale
            A = lambda v: sum(s * f.rmatvec(s * f.matvec(v)) for f, s in zip(facs, sw)) + mu * v
            x = torch.zeros_like(theta); rr = -g.clone(); p = rr.clone(); rs = rr @ rr
            g0 = rs.sqrt()
            for _ in range(cg_iters):
                Ap = A(p); alpha = rs / (p @ Ap)
                x = x + alpha * p; rr = rr - alpha * Ap
                rs_new = rr @ rr
                if rs_new.sqrt() < cg_tol * g0:
                    break
                p = rr + (rs_new / rs) * p; rs = rs_new
            L_new = total(theta + x)
            if L_new < L:
                theta, L, lam, accepted = theta + x, L_new, max(lam / 3.0, 1e-12), True
                break
            lam *= 4.0
        rec = {'iter': it, 'loss': L, 'lam': lam, 'accepted': accepted, 'elapsed_s': time.time() - t0}
        log.append(rec)
        if callback is not None:
            callback(it, theta.view(shape), rec)
        if not accepted:
            break
    return theta.view(shape), log
