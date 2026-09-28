"""Ambient CS fields and independent coordinate/finite-difference checks.

Local charts are used only for verification, not as the learned representation.
"""
from __future__ import annotations

import copy
import torch
from torch import nn

from extra.campaign_unit import Net, make_field
from src.loss.ambient_loss import compute_J_and_DJ, compute_P_and_DP, ambient_nijenhuis

DTYPE = torch.float64


class Correction(nn.Module):
    """A local ACS deformation; no orthogonality or analytic target is imposed."""
    def __init__(self, parent, base, geo, hidden):
        super().__init__()
        self.parent = parent.requires_grad_(False)
        self.head = Net(geo.N, geo.N, 1, hidden, .1).to(DTYPE)
        nn.init.zeros_(self.head.net[-1].weight)
        nn.init.zeros_(self.head.net[-1].bias)
        self.base, self.geo = base, geo
        self.stages = parent.stages + 1

    def forward(self, x):
        raw = self.parent(x)
        P, J = self.geo.manifold.projector(x), self.base(x)
        for stage in range(self.parent.stages):
            B = P @ raw[:, stage] @ P
            J = torch.matrix_exp(B) @ J @ torch.matrix_exp(-B)
        A = P @ self.head(x)[:, 0] @ P
        C = .5 * (A + J @ A @ J)
        return torch.cat((raw, C[:, None]), dim=1)


def architecture(ck):
    if ck.get('model_type') == 'precision':
        return ck['architecture']
    if ck.get('model_type') != 'mlp':
        raise ValueError('Expected an MLP or precision checkpoint')
    return {'kind': 'mlp', 'hidden': ck['hidden'], 'stages': ck['stages']}


def construct(arch, geo, base):
    if arch['kind'] == 'mlp':
        return Net(geo.N, geo.N, arch['stages'], arch['hidden'], .1).to(DTYPE)
    if arch['kind'] == 'correction':
        return Correction(construct(arch['parent'], geo, base), base, geo, arch['hidden'])
    raise ValueError('Unknown precision architecture')


def load_model(ck, geo, base):
    if any(ck.get(k) != v for k, v in {'case': 's4xs1xs1', 'task': 'cs', 'kind': 'ambient'}.items()):
        raise ValueError('Expected an ambient S4xS1xS1 CS checkpoint')
    arch = architecture(ck)
    with torch.random.fork_rng():
        model = construct(arch, geo, base)
    model.load_state_dict(ck['model'])
    return model, arch


def payload(model, arch, evidence, source):
    return {'case': 's4xs1xs1', 'task': 'cs', 'kind': 'ambient',
            'model_type': 'precision', 'architecture': copy.deepcopy(arch),
            'model': copy.deepcopy(model.state_dict()), 'stages': model.stages,
            'audited_base': evidence, 'source': source}


def tensor(geo, field, x):
    J, DJ = compute_J_and_DJ(field, x)
    P, DP = compute_P_and_DP(geo.manifold, x)
    return ambient_nijenhuis(J, DJ, P, DP)


def local_basis(point, dims):
    """Block-diagonal orthonormal basis, fixed while differentiating a chart."""
    E = point.new_zeros((len(point), sum(dims)))
    a = b = 0
    for d in dims:
        p = point[a:a+d+1]
        _, _, Vh = torch.linalg.svd(p[None], full_matrices=True)
        E[a:a+d+1, b:b+d] = Vh[1:].T
        a += d + 1
        b += d
    return E


def local_chart(field, point, dims, E):
    """q -> J in the chart normalize(point + E q); analytic chart differential."""
    def evaluate(q):
        blocks, differentials = [], []
        a = b = 0
        for d in dims:
            basis = E[a:a+d+1, b:b+d]
            y = point[a:a+d+1] + basis @ q[b:b+d]
            length = y.norm()
            p = y / length
            D = (torch.eye(d+1, dtype=q.dtype) - p[:, None]*p[None, :]) @ basis / length
            blocks.append(p)
            differentials.append(D)
            a += d + 1
            b += d
        x, D = torch.cat(blocks), torch.block_diag(*differentials)
        return torch.linalg.solve(D.T @ D, D.T @ field(x[None])[0] @ D)
    return evaluate


def coordinate_nijenhuis(J, D):
    """Independent coordinate formula: D[a,j,k] = partial_k J[a,j]."""
    first = torch.einsum('ki,ajk->aij', J, D)
    exterior = D.transpose(1, 2) - D
    return first - first.transpose(1, 2) - torch.einsum('ak,kij->aij', J, exterior)


def crosscheck(geo, field, points, check=lambda: None, tolerance=1e-5):
    rows = []
    for point in points:
        check()
        E = local_basis(point, geo.dims)
        chart = local_chart(field, point, geo.dims, E)
        q = point.new_zeros(geo.d)
        J, D = chart(q), torch.func.jacfwd(chart)(q)
        Nc = coordinate_nijenhuis(J, D)
        Na = tensor(geo, field, point[None])[0]
        converted = torch.einsum('ap,aij,iq,jr->pqr', E, Na, E, E)
        errors = []
        for h in (1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6):
            check()
            columns = []
            with torch.no_grad():
                for direction in torch.eye(geo.d, dtype=DTYPE):
                    columns.append((chart(h*direction)-chart(-h*direction))/(2*h))
                Nfd = coordinate_nijenhuis(J.detach(), torch.stack(columns, dim=-1))
                errors.append({'h': h, 'absolute_error': (Nfd-Nc).norm().item()})
        error = (Nc-converted).norm().item()
        row = {'point': point.tolist(), 'nij': Nc.norm().item(),
               'coordinate_ad_error': error, 'fd_sweep': errors,
               'covariant_DJ_norm': D.norm().item(),
               'tangent_J_condition': torch.linalg.cond(J).item()}
        finite = all(torch.isfinite(t).all().item() for t in (J, D, Nc, converted))
        row['passed'] = bool(finite and error <= 1e-8 + 1e-7*Nc.norm().item()
                             and sum(e['absolute_error'] <= tolerance for e in errors) >= 2)
        rows.append(row)
    return {'passed': bool(rows) and all(r['passed'] for r in rows), 'points': rows,
            'absolute_fd_tolerance': tolerance}


def six_torus_control():
    """Nonconstant pullback of an integrable T6 structure; diagnostic only."""
    from types import SimpleNamespace
    from src.geometry.ambient import AmbientManifold
    geo = SimpleNamespace(N=12, d=6, dims=(1,)*6, manifold=AmbientManifold((1,)*6))
    J0 = torch.block_diag(*[torch.tensor([[0., -1.], [1., 0.]], dtype=DTYPE)]*3)
    def field(x):
        # Avoid in-place operations on differentiable values under vmap.
        E = torch.stack([torch.nn.functional.pad(
            torch.stack((-x[:, 2*i+1], x[:, 2*i]), dim=1), (2*i, 10-2*i)) for i in range(6)], dim=1)
        F = torch.eye(6, dtype=DTYPE).expand(len(x), 6, 6)
        direction = torch.zeros((6, 6), dtype=DTYPE)
        direction[0, 1] = .2
        F = F + x[:, 2, None, None]*direction
        # The shear is nilpotent, so its inverse is exactly I - shear.
        # This also avoids platform-dependent batched solve forward-AD kernels.
        inverse = torch.eye(6, dtype=DTYPE) - x[:, 2, None, None]*direction
        J = inverse @ J0 @ F
        return E.transpose(1, 2) @ J @ E
    return geo, field
