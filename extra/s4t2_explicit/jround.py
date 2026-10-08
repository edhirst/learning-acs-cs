"""The complex structure J_round = F_* J_Y on the ROUND S^4 x T^2 = S^4 x S^1 x S^1 in R^5 x R^2 x R^2, as a function of
the ambient point q (docs/s4t2_explicit_diffeo.md, section 12.7 / 13).

    J, DJ, res = J_round(q)        q: (B, 9) points of S^4 x T^2 (unit S^4, unit circles)
        J   (B, 9, 9)     tangential complex structure (J = P J P, J^2 = -P), the layout of src/loss/ambient_loss.py
        DJ  (B, 9, 9, 9)  DJ[b, a, r, c] = d_c J_ar along tangent directions (zero on normals), exact (nested AD)
        res (B,)          |F(F^-1(q)) - q|_inf of the inverse (about 1e-13)

Evaluation: p = F^-1(q) in a holomorphic chart of Y chosen automatically (tinverse.Finv), J = DF J0 DF^+ there, DJ
from D^2 F (nijenhuis_round.transport). J_round is invariant under the T^3 = S^1 x T^2 action
(Z -> e^{i theta} Z, circles rotated), because F is T^3-equivariant and T^3 acts holomorphically on Y.

Known limitation (campaign T14): the inverse can fail on exact central-fibre inputs outside the linear caps.
The caps, including both exact poles, use a closed inverse and have dedicated J/DJ checks.
"""
from __future__ import annotations

import math
import torch

from extra.s4t2_explicit import nijenhuis_round as NR, tinverse as TI


def linear_ball_coframe(q, key):
    """Explicit real coframe matrix in (Re Z1,Im Z1,Re Z2,Im Z2,phi2,phi3).

    The complete formula is programme section 3.4 (C9). No automatic differentiation,
    inverse-map evaluation, angle branch or division by a vanishing chart coordinate.
    The exact definition uses the convergent infinite height series; J=8 below is
    the same finite evaluation used by tcore.mu_chart, not a tail certificate.
    """
    if key not in ('c0', 'ci'):
        raise ValueError('linear-ball coframe requires c0 or ci')
    C, scale, sigma = (TI.T.CUSP0, TI.L.S0, 1) if key == 'c0' else (TI.T.CUSPINF, TI.L.SINF, -1)
    radius = (q[:, :4]**2).sum(-1) / scale
    if ((sigma*q[:, 4] >= 0) | (radius > TI.L.EPS1*(1+1e-12))).any():
        raise ValueError('point outside the selected linear cap')
    Z1, Z2 = torch.complex(q[:,0],q[:,1]), torch.complex(q[:,2],q[:,3])
    dZ1 = torch.tensor([1,1j,0,0,0,0], dtype=Z1.dtype, device=q.device)
    dZ2 = torch.tensor([0,0,1,1j,0,0], dtype=Z1.dtype, device=q.device)
    phase = -1j if key == 'c0' else 1
    xs = ((Z1+phase*Z2)/math.sqrt(2*scale)).conj()/math.sqrt(C.a0)
    ys = (Z1-phase*Z2)/math.sqrt(2*scale)
    dxs = ((dZ1+phase*dZ2)/math.sqrt(2*scale)).conj()/math.sqrt(C.a0)
    dys = (dZ1-phase*dZ2)/math.sqrt(2*scale)
    t = xs*ys
    dt = ys[:,None]*dxs+xs[:,None]*dys
    X, Y = C.a0*TI.T.abs2(xs), TI.T.abs2(ys)
    dX = 2*C.a0*(xs.conj()[:,None]*dxs).real
    dY = 2*(ys.conj()[:,None]*dys).real
    tt = TI.T.abs2(t)
    f = 2*tt/(1-tt)
    df = 4*(t.conj()[:,None]*dt).real/(1-tt[:,None])**2
    num, den = torch.zeros_like(X), torch.zeros_like(X)
    dnum, dden = torch.zeros_like(dX), torch.zeros_like(dX)
    for j in range(-8,9):
        a, b, jj = j*(j+1)//2, j*(j-1)//2, j*(j+1)
        e = torch.exp(-jj*f)
        term = e*X**a*Y**b
        dterm = -jj*term[:,None]*df
        if a:
            dterm = dterm+(a*e*X**(a-1)*Y**b)[:,None]*dX
        if b:
            dterm = dterm+(b*e*X**a*Y**(b-1))[:,None]*dY
        num, den = num+j*term, den+term
        dnum, dden = dnum+j*dterm, dden+dterm
    H = num/den
    dH = (dnum-H[:,None]*dden)/den[:,None]
    dW = torch.tensor([0,0,0,0,1,-1j*sigma], dtype=Z1.dtype, device=q.device)
    # T^2 acts holomorphically, so the tensor in these coordinates is independent
    # of phi. Set W=phi2-i*sigma*phi3 to zero at the evaluation point, retaining dW.
    dA = t[:,None]*dW+H[:,None]*df+f[:,None]*dH
    eta1 = dxs-xs[:,None]*dA
    eta2 = dys+ys[:,None]*dA
    eta3 = sigma*(dW.conj()-t[:,None]**2*dW+2*H[:,None]*dt+2*t[:,None]*dH)
    return torch.stack([eta1.real,eta1.imag,eta2.real,eta2.imag,eta3.real,eta3.imag],-2)


def J_linear_ball(q, key):
    """Explicit coframe expression for J throughout either pole cap (value only).

    Symbolically J6=adj(M) J0 M/det(M), with
    det(M)=-sigma*(1-|t|^4)/(scale^2*a0). This numerical evaluation uses a 6x6
    inverse. The returned ambient 9x9 tensor annihilates all three normal directions.
    The existing global J_round path is kept as an independent comparison.
    """
    M = linear_ball_coframe(q,key)
    J6 = torch.linalg.inv(M) @ NR.J0.to(q) @ M
    E = torch.zeros(len(q),9,6,dtype=q.dtype,device=q.device)
    E[:,:4,:4] = torch.eye(4,dtype=q.dtype,device=q.device)
    E[:,4,:4] = -q[:,:4]/q[:,4,None]
    E[:,5,4], E[:,6,4] = -2*math.pi*q[:,6], 2*math.pi*q[:,5]
    E[:,7,5], E[:,8,5] = -2*math.pi*q[:,8], 2*math.pi*q[:,7]
    R = torch.zeros(len(q),6,9,dtype=q.dtype,device=q.device)
    R[:,:4,:5] = NR.S4T2.projector(q)[:,:4,:5]
    R[:,4,5], R[:,4,6] = -q[:,6]/(2*math.pi), q[:,5]/(2*math.pi)
    R[:,5,7], R[:,5,8] = -q[:,8]/(2*math.pi), q[:,7]/(2*math.pi)
    return E @ J6 @ R


def J_round(q, with_DJ=True):
    sol = TI.Finv(q)
    B = q.shape[0]
    J = torch.zeros(B, 9, 9)
    DJ = torch.zeros(B, 9, 9, 9) if with_DJ else None
    res = torch.zeros(B)
    for key, (idx, p, r) in sol.items():
        f = TI.CHART_F[key]
        if with_DJ:
            qq, A, H = NR.DF_and_D2F(f, p)
            Jk, DJk = NR.transport(A, H, qq)
            DJ[idx] = DJk
        else:
            A = torch.stack([torch.func.jvp(f, (p,), (NR.EYE6[k].expand_as(p),))[1] for k in range(6)], -1)
            Jk = A @ NR.J0 @ torch.linalg.pinv(A)
        J[idx] = Jk
        res[idx] = r
    return J, DJ, res


def elementary_cap_map(u, key):
    """G in a holomorphic cusp chart for the NEW parametrization (programme 9.3).

    u=(Re Z1,Im Z1,Re Z2,Im Z2,phi2,phi3), with |Z_old|^2<=.01.
    This finite elementary map uses the torus gauge C=phi on the caps. It is
    not the original J_round parametrization and contains no height series.
    """
    if key not in ('c0','ci'):
        raise ValueError('elementary cap requires c0 or ci')
    sigma = 1 if key == 'c0' else -1
    r = (u[...,:4]**2).sum(-1)
    X = -sigma*torch.sqrt(1-r)
    Z = torch.stack([torch.complex(u[...,0],u[...,1]),torch.complex(u[...,2],u[...,3])],-1)
    xs,ys = TI.algebraic_cap_slice(Z,X,infinity=key=='ci')
    t = xs*ys
    W = torch.complex(u[...,4],-sigma*u[...,5])
    x,y = xs*torch.exp(-t*W),ys*torch.exp(t*W)
    w = sigma*(W.conj()-W*t*t)
    return torch.stack([x.real,x.imag,y.real,y.imag,w.real,w.imag],-1)


def elementary_cap_coframe(q, key):
    """Finite elementary real coframe for the new tensor, including the pole.

    Analytic first derivatives, no AD, series, quadrature, or implicit root.
    Torus equivariance permits evaluation at phi=0 while retaining dphi.
    """
    if key not in ('c0','ci'):
        raise ValueError('elementary cap requires c0 or ci')
    C,scale,sigma = (TI.T.CUSP0,TI.L.S0,1) if key=='c0' else (TI.T.CUSPINF,TI.L.SINF,-1)
    r = (q[:,:4]**2).sum(-1)
    s = torch.sqrt(1-r*r)
    old_r = r/(1+s)
    if ((sigma*q[:,4]>=0)|(old_r > .01*(1+1e-12))).any():
        raise ValueError('point outside the selected elementary cap')
    gain = torch.rsqrt(1+s)
    dgain = r/(2*s*(1+s)**1.5)
    D = gain[:,None,None]*torch.eye(4,dtype=q.dtype,device=q.device)
    D = D+2*dgain[:,None,None]*q[:,:4,None]*q[:,None,:4]
    Z1 = torch.complex(q[:,0],q[:,1])*gain
    Z2 = torch.complex(q[:,2],q[:,3])*gain
    dZ1 = torch.cat([torch.complex(D[:,0,:],D[:,1,:]),torch.zeros(len(q),2,dtype=Z1.dtype,device=q.device)],-1)
    dZ2 = torch.cat([torch.complex(D[:,2,:],D[:,3,:]),torch.zeros(len(q),2,dtype=Z1.dtype,device=q.device)],-1)
    phase = -1j if key=='c0' else 1
    xs = ((Z1+phase*Z2)/math.sqrt(2*scale)).conj()/math.sqrt(C.a0)
    ys = (Z1-phase*Z2)/math.sqrt(2*scale)
    dxs = ((dZ1+phase*dZ2)/math.sqrt(2*scale)).conj()/math.sqrt(C.a0)
    dys = (dZ1-phase*dZ2)/math.sqrt(2*scale)
    t = xs*ys
    dW = torch.tensor([0,0,0,0,1,-1j*sigma],dtype=Z1.dtype,device=q.device)
    eta1 = dxs-xs[:,None]*t[:,None]*dW
    eta2 = dys+ys[:,None]*t[:,None]*dW
    eta3 = sigma*(dW.conj()-t[:,None]**2*dW)
    return torch.stack([eta1.real,eta1.imag,eta2.real,eta2.imag,eta3.real,eta3.imag],-2)


def J_elementary_cap(q, key):
    """New cap tensor E adj(M) J0 M R/det(M); not a global evaluator.

    Matrix inversion below is finite linear algebra (the exact formula uses
    adjugate/determinant), not the inverse of an implicitly specified function.
    """
    M = elementary_cap_coframe(q,key)
    J6 = torch.linalg.inv(M)@NR.J0.to(q)@M
    E = torch.zeros(len(q),9,6,dtype=q.dtype,device=q.device)
    E[:,:4,:4] = torch.eye(4,dtype=q.dtype,device=q.device)
    E[:,4,:4] = -q[:,:4]/q[:,4,None]
    E[:,5,4],E[:,6,4] = -2*math.pi*q[:,6],2*math.pi*q[:,5]
    E[:,7,5],E[:,8,5] = -2*math.pi*q[:,8],2*math.pi*q[:,7]
    R = torch.zeros(len(q),6,9,dtype=q.dtype,device=q.device)
    R[:,:4,:5] = NR.S4T2.projector(q)[:,:4,:5]
    R[:,4,5],R[:,4,6] = -q[:,6]/(2*math.pi),q[:,5]/(2*math.pi)
    R[:,5,7],R[:,5,8] = -q[:,8]/(2*math.pi),q[:,7]/(2*math.pi)
    return E@J6@R


def nijenhuis(q, J, DJ):
    """Per-sample intrinsic round-metric |N| and the cancellation scale (sum of the bracket-term norms)."""
    Nn, scale, _ = NR.nij(q, J, DJ)
    return Nn, scale
