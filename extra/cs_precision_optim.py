"""Chunked least squares optimizers with transactional steps and fixed weights."""
from __future__ import annotations

import math
import torch


class Parameters:
    def __init__(self, model):
        self.model = model
        self.items = [(name, p) for name, p in model.named_parameters() if p.requires_grad]

    def vector(self):
        return torch.cat([p.detach().reshape(-1) for _, p in self.items])

    def dictionary(self, vector):
        result, offset = {}, 0
        for name, p in self.items:
            result[name] = vector[offset:offset+p.numel()].view_as(p)
            offset += p.numel()
        return result

    def assign(self, vector):
        with torch.no_grad():
            for (_, p), value in zip(self.items, self.dictionary(vector).values()):
                p.copy_(value)


def weights(scores, uniform_count, tail_weight, tail_fraction=.1):
    """Uniform MSE plus an explicit top-tail MSE; fixed within an inner solve."""
    if not 0 < uniform_count <= len(scores) or tail_weight < 0:
        raise ValueError('Invalid collocation weighting')
    w = torch.zeros_like(scores)
    w[:uniform_count] = 1 / uniform_count
    k = max(1, math.ceil(tail_fraction*len(scores)))
    w[scores.topk(k).indices] += tail_weight/k
    return w


def loss_value(residual, vector, points, w, chunk, check):
    total = 0.
    for x, weight in zip(points.split(chunk), w.split(chunk)):
        check()
        with torch.no_grad():
            r = residual(vector, x)
            value = (r.square().sum(dim=1)*weight).sum()
        if not torch.isfinite(value):
            raise FloatingPointError('Nonfinite least squares objective')
        total += value.item()
    return total


def cg(matvec, rhs, iterations, tolerance, check):
    x = torch.zeros_like(rhs)
    r, p = rhs.clone(), rhs.clone()
    rr, initial = r.dot(r), r.norm()
    done = 0
    for _ in range(iterations):
        check()
        if r.norm() <= tolerance*initial or initial == 0:
            break
        Ap = matvec(p)
        denominator = p.dot(Ap)
        if not torch.isfinite(denominator) or denominator <= 0:
            raise FloatingPointError('Damped GN system lost positive definiteness')
        alpha = rr/denominator
        x, r = x+alpha*p, r-alpha*Ap
        next_rr = r.dot(r)
        p, rr = r+(next_rr/rr)*p, next_rr
        done += 1
    return x, done


def gn_step(residual, vector, points, w, *, damping, chunk, cg_steps, check,
            cg_tolerance=.05, backtracks=8):
    """Matrix-free GN on signed components, not scalar residual norms.

    No mutation occurs until the caller accepts the returned vector. JVP/VJP
    graphs are rebuilt per chunk to bound memory. W and collocation are fixed.
    """
    gradient = torch.zeros_like(vector)
    before = 0.
    for x, weight in zip(points.split(chunk), w.split(chunk)):
        check()
        f = lambda v: residual(v, x)*weight.sqrt()[:, None]
        r, pullback = torch.func.vjp(f, vector)
        gradient += pullback(r)[0].detach()
        before += r.square().sum().item()
        del r, pullback
    if not torch.isfinite(gradient).all() or not math.isfinite(before):
        raise FloatingPointError('Nonfinite GN linearization')

    def matvec(direction):
        output = damping*direction
        for x, weight in zip(points.split(chunk), w.split(chunk)):
            check()
            f = lambda v: residual(v, x)*weight.sqrt()[:, None]
            _, derivative = torch.func.jvp(f, (vector,), (direction,))
            _, pullback = torch.func.vjp(f, vector)
            output = output + pullback(derivative)[0].detach()
            del derivative, pullback
        return output

    delta, iterations = cg(matvec, -gradient, cg_steps, cg_tolerance, check)
    accepted, ratio, after, scale = False, 0., before, 0.
    for backtrack in range(backtracks):
        scale = 2.**(-backtrack)
        trial = vector + scale*delta
        try:
            after = loss_value(residual, trial, points, w, chunk, check)
        except FloatingPointError:
            after = float('inf')  # reject trial, not the last valid iterate
        predicted = -2*scale*gradient.dot(delta).item()
        if after <= before - 1e-4*predicted and after < before:
            # Actual/predicted *linear* decrease; an operational damping rule.
            ratio = (before-after)/max(predicted, 1e-300)
            accepted = True
            break
    new_damping = max(1e-9, damping*.5) if accepted and ratio > .25 else min(1e9, damping*4)
    return (trial.detach() if accepted else vector.clone()), {
        'before': before, 'after': after if accepted else before,
        'accepted': accepted, 'cg_iterations': iterations, 'step_scale': scale if accepted else 0.,
        'gradient_norm': gradient.norm().item(), 'damping': damping, 'next_damping': new_damping}


def backward_objective(residual, params, points, w, chunk, check):
    """Backpropagate a full, fixed objective through functional model parameters."""
    # A differentiable vector connects each functional call to model parameters.
    total = 0.
    for x, weight in zip(points.split(chunk), w.split(chunk)):
        check()
        vector = torch.cat([p.reshape(-1) for _, p in params.items])
        r = residual(vector, x)
        loss = (r.square().sum(dim=1)*weight).sum()
        if not torch.isfinite(loss):
            raise FloatingPointError('Nonfinite optimization objective')
        loss.backward()
        total += loss.detach().item()
    if any(p.grad is None or not torch.isfinite(p.grad).all() for _, p in params.items):
        raise FloatingPointError('Missing/nonfinite optimization gradients')
    return total
