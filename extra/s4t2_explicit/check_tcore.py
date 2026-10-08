"""tcore.py (torch) vs the verified numpy modules, on random points: every shared function to round-off."""
from __future__ import annotations
import math, os, sys
import numpy as np
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from extra.s4t2_explicit import tcore as T
from extra.s4t2_explicit import model as Mo, glue as G, hopf as H, curve as Cu, lam as La
from extra.s4t2_explicit.cusp import CUSP0, CUSPINF, bump_step

torch.set_default_dtype(torch.float64)
rng = np.random.default_rng(11)
n = 2000
tt = lambda a: torch.as_tensor(a)
rel = lambda a, b: float(np.max(np.abs(np.asarray(a) - np.asarray(b)) / (1 + np.abs(np.asarray(b)))))
worst = []


def rep(name, err, tol=1e-12):
    worst.append((name, err))
    print(f'{name:<44s} {err:.1e}  {"ok" if err < tol else "FAIL"}')


u = rng.uniform(-0.5, 1.5, n)
rep('bump_step', rel(T.bump_step(tt(u)).numpy(), bump_step(u)))
r = np.exp(rng.uniform(np.log(0.05), np.log(20), n)); th = rng.uniform(-np.pi, np.pi, n)
t = r * np.exp(1j * th); logt = np.log(r) + 1j * th
x = rng.uniform(-1, 2, (n, 4))
L = Mo.lattice(t)
z, w = (x * L[:, 0, :]).sum(-1), (x * L[:, 1, :]).sum(-1)
rep('lattice_coords', rel(T.lattice_coords(tt(z), tt(w), tt(t), tt(logt)).numpy(), x))
for C, TC, nm in ((CUSP0, T.CUSP0, 'c0'), (CUSPINF, T.CUSPINF, 'ci')):
    s = 0.6 * np.sqrt(rng.uniform(0, 1, n)) * np.exp(1j * rng.uniform(-np.pi, np.pi, n)) + 1e-12
    zeta = np.exp(rng.uniform(-4, 4, n) + 1j * rng.uniform(-np.pi, np.pi, n))
    rep(f'{nm}: P, multiplier', max(rel(TC.P(tt(s)).numpy(), C.P(s)), rel(TC.multiplier(tt(s)).numpy(), C.multiplier(s))))
    rep(f'{nm}: mu_torus', rel(TC.mu_torus(tt(zeta), tt(s)).numpy(), C.mu_torus(zeta, s)))
    rep(f'{nm}: E_torus', rel(TC.E_torus(tt(zeta), tt(s)).numpy(), C.E_torus(zeta, s)))
    xx = np.sqrt(rng.uniform(0, 1, n)) * np.exp(1j * rng.uniform(-np.pi, np.pi, n))     # |x|, |y| <= 1: the chart's
    yy = np.sqrt(rng.uniform(0, 1, n)) * np.exp(1j * rng.uniform(-np.pi, np.pi, n))     # range (|mu| <= 0.6)
    xx[: n // 10] = 0.0                                                                 # central fibre t = 0 (D_0)
    yy[n // 10: n // 5] = 0.0                                                           # central fibre t = 0 (D_1)
    rep(f'{nm}: mu_chart', rel(TC.mu_chart(tt(xx), tt(yy)).numpy(), C.mu_chart(xx, yy, 0)))
    m = (np.abs(C.mu_chart(xx, yy, 0)) < 0.6) & (np.abs(xx * yy) > 0)
    rep(f'{nm}: E_chart vs E_torus(y, xy)', rel(TC.E_chart(tt(xx[m]), tt(yy[m])).numpy(), C.E_torus(yy[m], xx[m] * yy[m])))
    m0 = np.abs(xx * yy) == 0
    print(f'      E_chart on the central fibre finite: {bool(np.isfinite(TC.E_chart(tt(xx[m0]), tt(yy[m0])).numpy()).all())}')
    fn = H.a1_c0 if nm == 'c0' else H.a2_ci
    tfn = T.a1_c0 if nm == 'c0' else T.a2_ci
    rep(f'{nm}: a1 / a2 (cusp)', rel(np.exp(2j * np.pi * tfn(tt(zeta), tt(s)).numpy()), np.exp(2j * np.pi * fn(zeta, s))))
    hf, thf = (H.h_B_c0, T.h_B_c0) if nm == 'c0' else (H.h_B_ci, T.h_B_ci)
    rep(f'{nm}: h_B', rel(thf(tt(zeta), tt(s)).numpy(), hf(zeta, s)))
r = np.exp(rng.uniform(np.log(0.5), np.log(2.0), n)); t = r * np.exp(1j * rng.uniform(-np.pi, np.pi, n))
x4 = rng.uniform(0, 1, n); x1 = rng.uniform(0, 1, n); logt = np.log(t)
rep('xi4', rel(np.exp(2j * np.pi * T.xi4(tt(t)).numpy()), np.exp(2j * np.pi * H.xi4(t))))
rep('h_B_mid', rel(T.h_B_mid(tt(t), tt(x4)).numpy(), H.h_B_mid(t, x4)))
r = np.exp(rng.uniform(np.log(0.05), np.log(0.7), n)); t = r * np.exp(1j * rng.uniform(-np.pi, np.pi, n)); logt = np.log(t)
rep('E0_W', rel(T.E0_W(tt(t), tt(x1), tt(x4), tt(logt)).numpy(), G.E0_W(t, x1, x4, logt)))
rep('Einf_W', rel(T.Einf_W(tt(1 / t), tt(x1), tt(x4), tt(-logt)).numpy(), G.Einf_W(1 / t, x1, x4, -logt)))
uu = rng.normal(size=(n, 4)); uu /= np.linalg.norm(uu, axis=1, keepdims=True); uc = uu[:, :2] + 1j * uu[:, 2:]
rep('moebius', rel(T.moebius(tt(uc)).numpy(), Cu.moebius(uc)))
ww = rng.normal(size=(n, 3)) * np.exp(rng.uniform(-3, 3, (n, 1)))
v1, X1 = T.P_final(tt(ww)); v0, X0 = La.P_final(ww)
rep('P_final', max(rel(v1.numpy(), v0), rel(X1.numpy(), X0)))
# Independent finite-height reference: sum 128 positive-tail seed terms in numpy.
# This is a numerical cross-check, not the geometric-tail proof in L13/N6.
def seed_numpy(v):
    reciprocal = v > 2
    w = np.where(reciprocal,1/np.where(reciprocal,v,1),v)
    transition = (w>1/16)&(w<.5)
    weight = bump_step(np.log(16*np.where(transition,w,1/16))/np.log(8))
    out = np.where(w>=.5,w/(1+w),w+weight*(w/(1+w)-w))
    return np.where(reciprocal,1-out,out)

from torch.func import jvp
for C in (T.FINITE_CUSP0,T.FINITE_CUSPINF):
    rr = np.exp(rng.uniform(-12,np.log(.5),3000))
    tv = rr*np.exp(1j*rng.uniform(-np.pi,np.pi,3000))
    qv = C.a0*rr**2*np.exp(-4*rr**2/(1-rr**2))
    yv = np.exp(.5*np.log(qv)*rng.uniform(-2,3,3000))
    xv = tv/yv
    av,y2v = qv/yv**2,yv**2
    ref = np.zeros(3000)
    for _ in range(128):
        ref += seed_numpy(av)-seed_numpy(y2v)
        av,y2v = qv*av,qv*y2v
    hc = C.mu_chart(tt(xv),tt(yv).to(torch.complex128))
    ht = C.mu_torus(tt(yv).to(torch.complex128),tt(tv))
    rep(f'finite {C.name}: 3000 chart/independent sums',rel(hc,ref))
    rep(f'finite {C.name}: chart/torus formula',rel(hc,ht))
    for shift in (-2,-1,1,2):
        hd = C.mu_torus(tt(yv).to(torch.complex128)*C.multiplier(tt(tv))**shift,tt(tv))
        rep(f'finite {C.name}: deck shift {shift}',rel(hd,ht+shift))
    yy = tt(yv).to(torch.complex128)
    def along(rho):
        return C.mu_torus(torch.exp(rho).to(torch.complex128),tt(tv))
    rho = yy.abs().log()
    deriv = jvp(along,(rho,),(torch.ones_like(rho),))[1]
    assert (deriv<0).all()
    # Exact axes, and derivatives through their zero coordinate. The finite
    # polynomial-in-q reference has the same first three jets, independently
    # of the implementation's batch-dependent geometric tail count.
    rad = torch.tensor([0.,.02,.05,1/16,.1,.5,1.,2.,16.,20.])
    zero = torch.zeros_like(rad).to(torch.complex128)
    axis = rad.sqrt().to(torch.complex128)
    rep(f'finite {C.name}: both exact axes',max(
        rel(C.mu_chart(axis/math.sqrt(C.a0),zero),T.finite_height_seed(rad)),
        rel(C.mu_chart(zero,axis),-T.finite_height_seed(rad))))
    def chart(z):
        xx,yy = torch.complex(z[:,0],z[:,1]),torch.complex(z[:,2],z[:,3])
        return C.mu_chart(xx,yy)
    def reference(z):
        xx,yy = torch.complex(z[:,0],z[:,1]),torch.complex(z[:,2],z[:,3])
        ts = T.abs2(xx*yy)
        m2 = C.a0*torch.exp(-4*ts/(1-ts))
        aa,bb,qq = m2*T.abs2(xx),T.abs2(yy),m2*ts
        out = torch.zeros_like(ts)
        for _ in range(8):
            out = out+T.finite_height_seed(aa)-T.finite_height_seed(bb)
            aa,bb = qq*aa,qq*bb
        return out
    z = torch.zeros(20,4)
    z[:10,0],z[10:,2] = rad.sqrt()/math.sqrt(C.a0),rad.sqrt()
    tangents = torch.tensor([[.23,-.17,.31,.11]]).expand_as(z)
    fa,fb = chart,reference
    for degree in range(1,4):
        pa,pb = fa,fb
        fa = lambda z,pa=pa:jvp(pa,(z,),(tangents,))[1]
        fb = lambda z,pb=pb:jvp(pb,(z,),(tangents,))[1]
        rep(f'finite {C.name}: axis/pole jet order {degree}',rel(fa(z),fb(z)))

bad = [nm for nm, e in worst if not e < 1e-12]
print('ALL PASS' if not bad else f'FAILURES: {bad}')
assert not bad,bad
