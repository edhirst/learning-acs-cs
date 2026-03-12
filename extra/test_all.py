"""
Comprehensive test suite for the two-patch ball-chart codebase.
Run with: pytest tests/test_all.py -v
"""

import torch
import numpy as np
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))

from src.geometry import EvenSphere, RADIAL_MIDPOINT
from src.models import TwoPatchAlmostComplexNet, PatchSubNet
from src.loss import NijenhuisTensorLoss, OverlapLoss, TotalLoss


class TestSphereGeometry:
    """Tests for ball-patch sphere geometry."""

    def test_sphere_creation(self):
        for dim in [2, 4, 6, 8]:
            sphere = EvenSphere(dimension=dim)
            assert sphere.dim == dim
            assert sphere.complex_dim == dim // 2
            assert sphere.ambient_dim == dim + 1

    def test_invalid_dimension(self):
        try:
            EvenSphere(dimension=3)
            assert False, "Should raise ValueError for odd dimension"
        except ValueError:
            pass

    def test_sample_points_patch1(self):
        sphere = EvenSphere(dimension=4)
        pts = sphere.sample_points_patch1(100)
        assert pts.shape == (100, 4)
        # All points must lie strictly inside the unit ball
        norms = torch.norm(pts, dim=1)
        assert (norms < 1.0 - 1e-6).all(), "Points must be strictly inside B(0,1)"

    def test_radial_midpoint(self):
        # RADIAL_MIDPOINT = √2 - 1
        assert abs(RADIAL_MIDPOINT - (2 ** 0.5 - 1)) < 1e-6

    def test_patch_change_involution(self):
        """patch_change_coords is self-inverse: φ(φ(x)) = x."""
        sphere = EvenSphere(dimension=4)
        x = sphere.sample_points_patch1(50)
        xx = EvenSphere.patch_change_coords(x)
        xxx = EvenSphere.patch_change_coords(xx)
        err = torch.norm(xxx - x, dim=1).max()
        assert err < 1e-5, f"Patch change not involutive: max err {err}"

    def test_patch_change_jacobian_shape(self):
        sphere = EvenSphere(dimension=4)
        x = sphere.sample_points_patch1(10)
        F = EvenSphere.patch_change_jacobian(x)
        assert F.shape == (10, 4, 4)

    def test_patch_change_J_preserves_structure(self):
        """J' = F J F⁻¹ still satisfies (J')² = −I if J² = −I."""
        dim = 4
        # Use fixed-radius equatorial points to avoid near-singular Jacobians (float32)
        torch.manual_seed(0)
        dirs = torch.randn(20, dim)
        dirs = dirs / dirs.norm(dim=1, keepdim=True)
        x = dirs * 0.4

        # Build J₀
        J0 = torch.zeros(dim, dim)
        for k in range(dim // 2):
            J0[2*k, 2*k+1] = -1.0
            J0[2*k+1, 2*k] = 1.0
        J = J0.unsqueeze(0).expand(20, -1, -1).clone()

        Jp = EvenSphere.patch_change_J(x, J)
        I  = torch.eye(dim).unsqueeze(0).expand(20, -1, -1)
        err = torch.norm(torch.bmm(Jp, Jp) + I, dim=(1, 2)).max()
        assert err < 1e-4, f"patch_change_J broke J²=−I: {err}"

    def test_ball_to_sphere(self):
        """ball_to_sphere outputs unit-norm ambient vectors."""
        sphere = EvenSphere(dimension=4)
        x = sphere.sample_points_patch1(30)
        X = EvenSphere.ball_to_sphere(x, from_north=True)
        assert X.shape == (30, 5)
        norms = torch.norm(X, dim=1)
        assert torch.allclose(norms, torch.ones(30), atol=1e-5)


class TestAlmostComplexNet:
    """Tests for the two-patch network."""

    def test_network_creation(self):
        net = TwoPatchAlmostComplexNet(
            manifold_dim=4,
            hidden_layers=[64, 128, 64],
            activation='tanh',
        )
        assert net.manifold_dim == 4
        assert net.complex_dim == 2
        assert net.num_dof == 16
        assert isinstance(net.patch1, PatchSubNet)
        assert isinstance(net.patch2, PatchSubNet)

    def test_forward_pass_shapes(self):
        net = TwoPatchAlmostComplexNet(
            manifold_dim=4,
            hidden_layers=[64, 128, 64],
            activation='tanh',
        )
        sphere = EvenSphere(dimension=4)
        x1 = sphere.sample_points_patch1(8)
        J1, J2, x1_out, x2 = net(x1)
        assert J1.shape == (8, 4, 4)
        assert J2.shape == (8, 4, 4)
        assert x1_out.shape == (8, 4)
        assert x2.shape == (8, 4)

    def test_structure_constraint(self):
        """J² = −I must hold exactly for both patches."""
        net = TwoPatchAlmostComplexNet(
            manifold_dim=4,
            hidden_layers=[64, 128, 64],
            activation='tanh',
        )
        sphere = EvenSphere(dimension=4)
        x1 = sphere.sample_points_patch1(16)
        J1, J2, x1, x2 = net(x1)

        I = torch.eye(4).unsqueeze(0).expand(16, -1, -1)
        err1 = torch.norm(torch.bmm(J1, J1) + I, dim=(1, 2)).mean()
        err2 = torch.norm(torch.bmm(J2, J2) + I, dim=(1, 2)).mean()
        assert err1 < 1e-5, f"J1²+I non-zero: {err1}"
        assert err2 < 1e-5, f"J2²+I non-zero: {err2}"

    def test_fourier_features(self):
        net = TwoPatchAlmostComplexNet(
            manifold_dim=4,
            hidden_layers=[64],
            activation='tanh',
            use_fourier_features=True,
            num_fourier_features=32,
        )
        sphere = EvenSphere(dimension=4)
        x1 = sphere.sample_points_patch1(8)
        J1, J2, _, _ = net(x1)
        assert J1.shape == (8, 4, 4)

    def test_different_activations(self):
        for act in ['tanh', 'relu', 'gelu', 'sin']:
            net = TwoPatchAlmostComplexNet(
                manifold_dim=2,
                hidden_layers=[32, 32],
                activation=act,
            )
            sphere = EvenSphere(dimension=2)
            x1 = sphere.sample_points_patch1(4)
            J1, J2, _, _ = net(x1)
            assert J1.shape == (4, 2, 2)


class TestLossFunctions:
    """Tests for Nijenhuis + overlap loss functions."""

    def _make_net_and_data(self, dim=4, batch=8):
        net = TwoPatchAlmostComplexNet(
            manifold_dim=dim,
            hidden_layers=[32, 64, 32],
            activation='tanh',
        )
        sphere = EvenSphere(dimension=dim)
        x1 = sphere.sample_points_patch1(batch)
        return net, x1

    def test_nijenhuis_loss_shape(self):
        net, x1 = self._make_net_and_data()
        loss_fn = NijenhuisTensorLoss(manifold_dim=4)
        loss, info = loss_fn(x1, net.patch1)
        assert loss.shape == ()
        assert 'nij_norm' in info
        assert loss.item() >= 0

    def test_overlap_loss_zero_for_consistent_J(self):
        """If J2 = F J1 F^{-1} exactly, overlap loss should be zero."""
        dim = 4
        batch = 16
        # Place all samples near the equator (rho ≈ 0.4) so the mask passes
        # every sample for the default r_min = 0.15.
        torch.manual_seed(0)
        dirs = torch.randn(batch, dim)
        dirs = dirs / dirs.norm(dim=1, keepdim=True)
        x1 = dirs * 0.4   # all radii exactly 0.4
        x2 = EvenSphere.patch_change_coords(x1)

        # Build J₀ exactly satisfying J² = −I
        J0 = torch.zeros(dim, dim)
        for k in range(dim // 2):
            J0[2*k, 2*k+1] = -1.0
            J0[2*k+1, 2*k] = 1.0
        # J1_in_2 = F J0 F^{-1} — pass this as J2 to get zero overlap
        J1 = J0.unsqueeze(0).expand(batch, -1, -1)
        J2 = EvenSphere.patch_change_J(x1, J1)

        loss_fn = OverlapLoss()
        loss, info = loss_fn(J1, J2, x1, x2)
        # One direction is exact; the other (patch-2→patch-1) should also be ~0
        assert loss.item() < 1e-4, f"Overlap loss should be zero: {loss.item()}"

    def test_total_loss(self):
        net, x1 = self._make_net_and_data()
        J1, J2, x1, x2 = net(x1)
        loss_fn = TotalLoss(manifold_dim=4, weights={'nijenhuis': 1.0, 'overlap': 1.0})
        loss, info = loss_fn(J1, J2, x1, x2, net)
        assert 'total_loss' in info
        assert 'nij1' in info
        assert 'nij2' in info
        assert 'overlap_norm' in info
        assert loss.item() >= 0

    def test_total_loss_backward(self):
        net, x1 = self._make_net_and_data()
        J1, J2, x1, x2 = net(x1)
        loss_fn = TotalLoss(manifold_dim=4, weights={'nijenhuis': 1.0, 'overlap': 1.0})
        loss, _ = loss_fn(J1, J2, x1, x2, net)
        loss.backward()
        # At least one parameter should have a non-None gradient
        grads = [p.grad for p in net.parameters() if p.grad is not None]
        assert len(grads) > 0, "No gradients computed"


class TestIntegration:
    """Integration tests: full training step."""

    def test_single_training_step(self):
        dim = 4
        sphere = EvenSphere(dimension=dim)
        net = TwoPatchAlmostComplexNet(
            manifold_dim=dim,
            hidden_layers=[32, 64, 32],
            activation='tanh',
        )
        loss_fn = TotalLoss(manifold_dim=dim, weights={'nijenhuis': 1.0, 'overlap': 1.0})
        optimizer = torch.optim.Adam(net.parameters(), lr=1e-3)

        x1 = sphere.sample_points_patch1(16)

        optimizer.zero_grad()
        J1, J2, x1, x2 = net(x1)
        loss, info = loss_fn(J1, J2, x1, x2, net)
        loss.backward()
        optimizer.step()

        assert not torch.isnan(loss), "Loss is NaN"
        assert not torch.isinf(loss), "Loss is Inf"

    def test_different_dimensions(self):
        for dim in [2, 4, 6]:
            sphere = EvenSphere(dimension=dim)
            net = TwoPatchAlmostComplexNet(
                manifold_dim=dim,
                hidden_layers=[32, 64, 32],
                activation='tanh',
            )
            x1 = sphere.sample_points_patch1(8)
            J1, J2, x1, x2 = net(x1)
            assert J1.shape == (8, dim, dim)
            assert J2.shape == (8, dim, dim)


def run_tests():
    """Run all tests manually (without pytest)."""
    print("Running tests...\n")

    test_classes = [
        TestSphereGeometry,
        TestAlmostComplexNet,
        TestLossFunctions,
        TestIntegration,
    ]

    total = passed = 0
    for cls in test_classes:
        print(f"\n{'='*60}")
        print(f"Testing: {cls.__name__}")
        print('='*60)
        instance = cls()
        methods = sorted(m for m in dir(instance) if m.startswith('test_'))
        for name in methods:
            total += 1
            try:
                getattr(instance, name)()
                print(f"  \u2713 {name}")
                passed += 1
            except Exception as exc:
                print(f"  \u2717 {name}: {exc}")

    print(f"\n{'='*60}")
    print(f"Results: {passed}/{total} tests passed")
    print('='*60)
    return passed == total


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
