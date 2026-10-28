"""T^2-invariant MLP conjugation on the round S^4 x T^2 (reduced model of
extra/s4t2_invariant.py): J = g2 g1 J_b g1^-1 g2^-1, g_s = exp(P B_s(x) P),
B_s from an MLP of x in S^4 only, J_b the unit's own exact T^2-invariant ACS
(degree-1 class from a random start, as in s4t2_invariant_run.py). Loss:
mean-square Nijenhuis norm with tail-aware sampling; Adam. Exact ACS by
construction; reports held-out mean/max and the adversarial sup of |N|.

    .venv/bin/python experiments/s4t2_invariant_conj.py --out outputs/paper_rerun/runs/s4t2invconj_s7 --seed 7
"""
from __future__ import annotations

import os
import sys
import json
import time
import argparse

import torch
import torch.nn as nn

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extra.s4t2_invariant import MAN, K_of, polar, nij_norm, adversarial_sup, acs_residual
from s4t2_invariant_run import acs_stage

DTYPE = torch.float64


class ConjNet(nn.Module):
    def __init__(self, hidden=(128, 128, 128), stages=2, scale=0.1):
        super().__init__()
        layers, prev = [], 5
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.GELU()]
            prev = h
        layers.append(nn.Linear(prev, stages * 49))
        self.net = nn.Sequential(*layers).to(DTYPE)
        with torch.no_grad():
            self.net[-1].weight.mul_(scale); self.net[-1].bias.normal_(0.0, scale)
        self.stages = stages

    def forward(self, p):
        return self.net(p[:, :5]).reshape(-1, self.stages, 7, 7)


def make_J(model, Jb):
    def J(p):
        P = MAN.projector(p)
        raw = model(p)
        Jx = Jb(p)
        for s in range(model.stages):
            B = P @ raw[:, s] @ P
            Jx = torch.matrix_exp(B) @ Jx @ torch.matrix_exp(-B)
        return Jx
    return J


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--steps', type=int, default=4000)
    ap.add_argument('--batch', type=int, default=128)
    ap.add_argument('--lr', type=float, default=1e-3)
    ap.add_argument('--init-from', default=None, help='run dir whose final.pt warm-starts the conjugator')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', 2)))
    t0 = time.time()
    logf = open(os.path.join(a.out, 'log.txt'), 'w')

    def log(msg):
        print(msg, flush=True); logf.write(msg + '\n'); logf.flush()
    fld, ladder = acs_stage(a.seed, log)
    Kb = K_of(fld)
    for p in fld.parameters():
        p.requires_grad_(False)
    Jb = polar(Kb)
    torch.manual_seed(a.seed)
    model = ConjNet()
    if a.init_from:
        model.load_state_dict(torch.load(os.path.join(a.init_from, 'final.pt')))
        log(f"warm start from {a.init_from}")
    J = make_J(model, Jb)
    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=a.steps, eta_min=1e-5)
    hist = []
    xv = MAN.sample(1000)
    for step in range(a.steps):
        with torch.no_grad():                                      # tail-aware batch
            pool = MAN.sample(2 * a.batch)
            sc = nij_norm(J, pool).detach()
        x = torch.cat([MAN.sample(a.batch // 2), pool[torch.topk(sc, a.batch // 2).indices]], 0)
        loss = (nij_norm(J, x) ** 2).mean()
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sch.step()
        if (step + 1) % 250 == 0 or step == 0:
            nv = torch.cat([nij_norm(J, xv[i:i + 250]).detach() for i in range(0, 1000, 250)])
            rec = {'step': step + 1, 'batch_loss': loss.item(), 'N_mean': nv.mean().item(), 'N_max': nv.max().item(),
                   'elapsed_min': (time.time() - t0) / 60}
            hist.append(rec)
            log(f"step {step+1:5d} | batch {loss.item():.3e} | held-out mean {rec['N_mean']:.3e} max {rec['N_max']:.3e} "
                f"| {rec['elapsed_min']:.1f} min")
    xs = MAN.sample(3000)
    nn_ = torch.cat([nij_norm(J, xs[i:i + 300]).detach() for i in range(0, 3000, 300)])
    adv = adversarial_sup(lambda p: nij_norm(J, p), 64, 60)
    with torch.no_grad():
        acs = acs_residual(J, xs[:500]).max().item()
    res = {'seed': a.seed, 'init_from': a.init_from, 'acs_ladder': ladder, 'hist': hist, 'N_mean': nn_.mean().item(), 'N_max': nn_.max().item(),
           'N_adv': adv, 'acs_residual_max': acs, 'steps': a.steps, 'elapsed_min': (time.time() - t0) / 60}
    log(f"RESULT: N mean {res['N_mean']:.3e} max {res['N_max']:.3e} adversarial {adv:.3e} (|J^2+P| {acs:.1e})")
    torch.save(model.state_dict(), os.path.join(a.out, 'final.pt'))
    json.dump(res, open(os.path.join(a.out, 'done.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
