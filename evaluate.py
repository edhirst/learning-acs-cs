"""
Utility script to evaluate a trained model checkpoint.
"""

import torch
import argparse

from src.geometry import EvenSphere, RADIAL_MIDPOINT
from src.geometry.sphere_product import SphereProduct
from src.models import TwoPatchAlmostComplexNet, FourPatchProductNet
from src.loss import TotalLoss
from src.loss.product_loss import ProductTotalLoss
from src.utils.chern import compute_first_chern_number


def _factored_mixing(J_dict: dict, dim_alpha: int) -> float:
    """
    Compute the mean factored-ACS mixing ratio across all four patches.

    For each patch, the mixing ratio is:

        mix(J) = (‖J_αβ‖_F² + ‖J_βα‖_F²) / ‖J‖_F²

    where J_αβ and J_βα are the off-diagonal blocks (α→β and β→α coupling).
    A value of 0 means J is block-diagonal (product structure); larger values
    indicate genuine coupling between the two sphere factors.

    Args:
        J_dict: {(i_α, i_β) → J tensor (batch, d, d)}.
        dim_alpha: m_α (size of the α block).

    Returns:
        Mean mixing ratio across all four patches (scalar float).
    """
    ratios = []
    for J in J_dict.values():
        ma = dim_alpha
        J_ab = J[:, :ma, ma:]
        J_ba = J[:, ma:, :ma]
        off_sq = (J_ab ** 2).sum(dim=(1, 2)) + (J_ba ** 2).sum(dim=(1, 2))
        tot_sq = (J  ** 2).sum(dim=(1, 2)).clamp(min=1e-12)
        ratios.append((off_sq / tot_sq).mean().item())
    return sum(ratios) / len(ratios)


def evaluate_model(
    checkpoint_path: str,
    num_samples: int = 1000,
):
    """
    Evaluate a trained model and print loss metrics.

    For sphere_product checkpoints, also prints the factored-ACS mixing ratio.

    Args:
        checkpoint_path: Path to model checkpoint.
        num_samples: Number of collocation points for evaluation.
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    config = checkpoint['config']

    print(f"Checkpoint epoch : {checkpoint['epoch']}")
    print(f"Best training loss: {checkpoint['best_loss']:.6e}")

    train_cfg = config.get('training', {})
    if 'precision' in train_cfg:
        precision = int(train_cfg['precision'])
    else:
        dtype_str = train_cfg.get('dtype', 'float32')
        precision = 64 if str(dtype_str).endswith('64') else 32
    dtype = torch.float64 if precision == 64 else torch.float32
    if dtype == torch.float64 and device.type == 'mps':
        print("Warning: MPS does not support float64; falling back to CPU.")
        device = torch.device('cpu')
    print(f"Using precision : float{precision}")

    manifold_type = config['geometry'].get('manifold_type', 'sphere')

    if manifold_type == 'sphere_product':
        dim_alpha = config['geometry']['real_dimension_alpha']
        dim_beta  = config['geometry']['real_dimension_beta']
        sphere_product = SphereProduct(dim_alpha, dim_beta, device=str(device))

        model = FourPatchProductNet(
            dim_alpha=dim_alpha,
            dim_beta=dim_beta,
            hidden_layers=config['network']['hidden_layers'],
            activation=config['network']['activation'],
            use_fourier_features=config['network']['use_fourier_features'],
            num_fourier_features=config['network']['num_fourier_features'],
            fourier_scale=config['network']['fourier_scale'],
            initialization=config['network']['initialization'],
            matrix_exp_terms=config['network'].get('matrix_exp_terms'),
            use_residual_blocks=config['network'].get('use_residual_blocks', False),
        ).to(device).to(dtype)
        model.load_state_dict(checkpoint['model_state_dict'])
        model.eval()

        loss_fn = ProductTotalLoss(
            dim_alpha=dim_alpha,
            dim_beta=dim_beta,
            weights=config['loss_weights'],
        ).to(device).to(dtype)

        xa, xb = sphere_product.sample_points(num_samples)
        xa, xb = xa.to(device=device, dtype=dtype), xb.to(device=device, dtype=dtype)
        J_dict, xa, xb, xa2, xb2 = model(xa, xb)

        osa = float(config['loss_weights'].get('overlap_size_alpha', 0.1))
        osb = float(config['loss_weights'].get('overlap_size_beta',  0.1))
        r_cut_a  = RADIAL_MIDPOINT + osa
        r_cut_b  = RADIAL_MIDPOINT + osb
        r_inner_a = (1.0 - r_cut_a) / (1.0 + r_cut_a)
        r_inner_b = (1.0 - r_cut_b) / (1.0 + r_cut_b)
        with torch.no_grad():
            prho_a = {0: xa.norm(dim=1),  1: xa2.norm(dim=1)}
            prho_b = {0: xb.norm(dim=1),  1: xb2.norm(dim=1)}
            n_nij = {
                (ia, ib): int(((prho_a[ia] <= r_cut_a) & (prho_b[ib] <= r_cut_b)).sum())
                for ia, ib in [(0, 0), (0, 1), (1, 0), (1, 1)]
            }
            ov_a    = (prho_a[0] >= r_inner_a) & (prho_a[0] <= r_cut_a)
            ov_b    = (prho_b[0] >= r_inner_b) & (prho_b[0] <= r_cut_b)
            n_ov_a    = int(ov_a.sum())
            n_ov_b    = int(ov_b.sum())
            n_ov_both = int((ov_a & ov_b).sum())

        _, info = loss_fn(J_dict, xa, xb, xa2, xb2, model, hard_cutoff=True, report_all=True)

        mixing = _factored_mixing(J_dict, dim_alpha)
        c1 = compute_first_chern_number(
            model,
            manifold_type='sphere_product',
            dim_alpha=dim_alpha,
            dim_beta=dim_beta,
            x_alpha=xa,
            x_beta=xb,
        )
        info['c1'] = c1

        d = dim_alpha + dim_beta
        print(f"\nEvaluation on {num_samples} samples (S^{dim_alpha} × S^{dim_beta}, d={d}):")
        for ia, ib in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            n = n_nij[(ia, ib)]
            print(f"  Nij {ia}{ib} region  : {n}/{num_samples} ({100*n/num_samples:.1f}%)")
        print(f"  Overlap α-band   : {n_ov_a}/{num_samples} ({100*n_ov_a/num_samples:.1f}%)  [α-flip pairs]")
        print(f"  Overlap β-band   : {n_ov_b}/{num_samples} ({100*n_ov_b/num_samples:.1f}%)  [β-flip pairs]")
        print(f"  Overlap αβ-band  : {n_ov_both}/{num_samples} ({100*n_ov_both/num_samples:.1f}%)  [both-flip pairs]")
        print(f"  Total loss    : {info['total_loss']:.6e}")
        for key in ['nij_00', 'nij_01', 'nij_10', 'nij_11']:
            print(f"  {key}         : {info.get(key, 0.0):.6e}")
        print(f"  Overlap       : {info['overlap_norm']:.6e}")
        print(f"  ACS mixing    : {mixing:.4f}  "
              f"(0 = block-diagonal / product structure, 1 = fully coupled)")
        print(f"  c₁           : {c1:.6f}")

    else:
        sphere = EvenSphere(dimension=config['geometry']['real_dimension'], device=str(device))

        model = TwoPatchAlmostComplexNet(
            manifold_dim=sphere.dim,
            hidden_layers=config['network']['hidden_layers'],
            activation=config['network']['activation'],
            use_fourier_features=config['network']['use_fourier_features'],
            num_fourier_features=config['network']['num_fourier_features'],
            fourier_scale=config['network']['fourier_scale'],
            initialization=config['network']['initialization'],
            matrix_exp_terms=config['network'].get('matrix_exp_terms'),
            use_residual_blocks=config['network'].get('use_residual_blocks', False),
        ).to(device).to(dtype)
        model.load_state_dict(checkpoint['model_state_dict'])
        model.eval()

        loss_fn = TotalLoss(manifold_dim=sphere.dim, weights=config['loss_weights']).to(device).to(dtype)

        x = sphere.sample_points_patch1(num_samples).to(device=device, dtype=dtype)
        J1, J2, x1, x2 = model(x)

        overlap_size = float(config['loss_weights'].get('overlap_size', 0.1))
        r_cut   = RADIAL_MIDPOINT + overlap_size
        r_inner = (1.0 - r_cut) / (1.0 + r_cut)
        with torch.no_grad():
            rho1 = x1.norm(dim=1)
            rho2 = x2.norm(dim=1)
            n_nij1 = int((rho1 <= r_cut).sum())
            n_nij2 = int((rho2 <= r_cut).sum())
            n_ov   = int(((rho1 >= r_inner) & (rho1 <= r_cut)).sum())

        _, info = loss_fn(J1, J2, x1, x2, model, hard_cutoff=True, report_all=True)
        c1 = compute_first_chern_number(
            model,
            manifold_type='sphere',
            manifold_dim=sphere.dim,
            x=x,
        )
        info['c1'] = c1

        print(f"\nEvaluation on {num_samples} samples (S^{sphere.dim}):")
        print(f"  Nij1 region  : {n_nij1}/{num_samples} ({100*n_nij1/num_samples:.1f}%)")
        print(f"  Nij2 region  : {n_nij2}/{num_samples} ({100*n_nij2/num_samples:.1f}%)")
        print(f"  Overlap band : {n_ov}/{num_samples} ({100*n_ov/num_samples:.1f}%)")
        print(f"  Total loss : {info['total_loss']:.6e}")
        print(f"  Nij1       : {info['nij1']:.6e}")
        print(f"  Nij2       : {info['nij2']:.6e}")
        print(f"  Overlap    : {info['overlap_norm']:.6e}")
        print(f"  c₁         : {c1:.6f}")

    return info


def main():
    parser = argparse.ArgumentParser(description='Evaluate trained model')
    parser.add_argument('checkpoint', type=str, help='Path to checkpoint')
    parser.add_argument('--samples', type=int, default=1000,
                        help='Number of samples for evaluation')
    args = parser.parse_args()

    evaluate_model(args.checkpoint, args.samples)


if __name__ == '__main__':
    main()
