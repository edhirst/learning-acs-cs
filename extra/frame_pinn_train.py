"""From-scratch frame-PINN: a network trained on the frame Nijenhuis residual
alone finds integrable complex structures on the (formerly failing) cases.

Setup is deliberately non-trivial:
  * the base structure J₀ᵣ = e^{A_r} J₀ e^{−A_r} is a FIXED RANDOM constant —
    generically NOT integrable (its loss is O(1−10)), so the network cannot
    succeed by outputting zero;
  * the network's output layer is initialised at scale O(1), so training
    starts from a random non-integrable FIELD;
  * the only training signal is the frame Nijenhuis residual (PINN loss).

The network must construct a field A(x) with e^{A(x)} J₀ᵣ e^{−A(x)}
integrable. Solutions exist (the integrable constants are in the orbit), but
the optimiser is told nothing about them.

Run:  python extra/frame_pinn_train.py [case]     (default s3xs3)
"""
from __future__ import annotations

import os
import sys
import time
import torch
import torch.nn as nn

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.geometry.frames import FRAMES
from src.loss.frame_loss import frame_nij_norm

torch.manual_seed(7)
DTYPE = torch.float64

ADAM_STEPS = 1500
BATCH = 256
LBFGS_ITERS = 150
HIDDEN = [128, 128, 128]


class FrameSubNet(nn.Module):
    """MLP: ambient point → generator A(x) ∈ ℝ^{d×d}; J = e^A J₀ᵣ e^{−A}."""

    def __init__(self, ambient_dim: int, d: int, J0r: torch.Tensor,
                 init_scale: float = 1.0):
        super().__init__()
        layers, prev = [], ambient_dim
        for h in HIDDEN:
            layers += [nn.Linear(prev, h), nn.GELU()]
            prev = h
        layers += [nn.Linear(prev, d * d)]
        self.net = nn.Sequential(*layers)
        with torch.no_grad():
            self.net[-1].weight.mul_(init_scale)
            self.net[-1].bias.normal_(0.0, init_scale)
        self.d = d
        self.register_buffer('J0r', J0r)

    def forward_A(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).reshape(-1, self.d, self.d)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        A = self.forward_A(x)
        return torch.matrix_exp(A) @ self.J0r @ torch.matrix_exp(-A)


def J_and_dJframe(model, x, E):
    """J (B,d,d) and frame-directional dJ (B,d,d,d) via forward AD + contraction."""
    def J_single(xs):
        return model(xs.unsqueeze(0)).squeeze(0)
    J = model(x)
    dJ_amb = torch.func.vmap(torch.func.jacfwd(J_single))(x)     # (B, d, d, N)
    dJ_frame = torch.einsum('beaN,bcN->beac', dJ_amb, E)         # deriv dir last
    return J, dJ_frame


def run(case: str):
    frame = FRAMES[case]()
    d, N = frame.dim, frame.ambient_dim
    c_const = None
    x_probe = frame.sample(200)
    c_probe = frame.structure_functions(x_probe)
    if (c_probe - c_probe.mean(dim=0, keepdim=True)).abs().max() < 1e-10:
        c_const = c_probe[:1]

    # canonical J0, then a fixed random conjugation to make it non-integrable
    J0 = torch.zeros(d, d, dtype=DTYPE)
    for k in range(d // 2):
        J0[2*k, 2*k+1] = -1.0
        J0[2*k+1, 2*k] = 1.0
    A_r = 0.6 * torch.randn(d, d, dtype=DTYPE)
    J0r = torch.matrix_exp(A_r) @ J0 @ torch.matrix_exp(-A_r)

    def c_of(x):
        if c_const is not None:
            return c_const.expand(len(x), d, d, d)
        return frame.structure_functions(x)

    # how bad is the base structure?
    xs = frame.sample(500)
    base_nij = frame_nij_norm(J0r.unsqueeze(0).expand(500, d, d),
                              torch.zeros(500, d, d, d, dtype=DTYPE), c_of(xs))
    print(f"[{case}] base J0r (random constant): frame ||N|| mean {base_nij.mean():.4f}"
          f"  (non-integrable start — outputting A=0 is NOT a solution)")

    model = FrameSubNet(N, d, J0r, init_scale=1.0).to(DTYPE)
    n_par = sum(p.numel() for p in model.parameters())

    x0 = frame.sample(500)
    with torch.no_grad():
        E0 = frame.frame_fields(x0)
    J_init, dJ_init = J_and_dJframe(model, x0, E0)
    init_loss = frame_nij_norm(J_init.detach(), dJ_init.detach(), c_of(x0)).mean()
    print(f"[{case}] network init field:          frame ||N|| mean {init_loss:.4f}"
          f"   ({n_par} params)")

    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=ADAM_STEPS, eta_min=1e-5)
    t0 = time.time()
    for step in range(ADAM_STEPS):
        x = frame.sample(BATCH)
        E = frame.frame_fields(x)
        J, dJ = J_and_dJframe(model, x, E)
        loss = (frame_nij_norm(J, dJ, c_of(x)) ** 2).mean()
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sched.step()
        if step % 250 == 0 or step == ADAM_STEPS - 1:
            print(f"  adam step {step:5d} | mean ||N||^2 {loss.item():.6e}")

    # LBFGS polish on a fixed batch
    xf = frame.sample(1024)
    Ef = frame.frame_fields(xf)
    cf = c_of(xf)
    lb = torch.optim.LBFGS(model.parameters(), lr=1.0, max_iter=LBFGS_ITERS,
                           history_size=100, line_search_fn='strong_wolfe')

    def closure():
        lb.zero_grad()
        J, dJ = J_and_dJframe(model, xf, Ef)
        l = (frame_nij_norm(J, dJ, cf) ** 2).mean()
        l.backward()
        return l
    lb.step(closure)

    # held-out evaluation
    xt = frame.sample(4000)
    Et = frame.frame_fields(xt)
    with torch.no_grad():
        pass
    J, dJ = J_and_dJframe(model, xt, Et)
    nij = frame_nij_norm(J.detach(), dJ.detach(), c_of(xt))
    print(f"[{case}] TRAINED, held-out 4000 pts:  frame ||N|| mean {nij.mean():.6e} "
          f"(max {nij.max():.2e})   [{time.time()-t0:.0f}s]")
    print(f"[{case}] path: random non-integrable field ({init_loss:.2f}) "
          f"→ integrable ({nij.mean():.2e}) by PINN residual alone")
    return model


if __name__ == '__main__':
    case = sys.argv[1] if len(sys.argv) > 1 else 's3xs3'
    run(case)
