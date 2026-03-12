"""
Utility script to evaluate a trained model checkpoint.
"""

import torch
import argparse

from src.geometry import EvenSphere
from src.geometry.sphere_product import SphereProduct
from src.models import TwoPatchAlmostComplexNet, FourPatchProductNet
from src.loss import TotalLoss
from src.loss.product_loss import ProductTotalLoss


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


def evaluate_model(checkpoint_path: str, num_samples: int = 1000):
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
        ).to(device)
        model.load_state_dict(checkpoint['model_state_dict'])
        model.eval()

        loss_fn = ProductTotalLoss(
            dim_alpha=dim_alpha,
            dim_beta=dim_beta,
            weights=config['loss_weights'],
        ).to(device)

        xa, xb = sphere_product.sample_points(num_samples)
        xa, xb = xa.to(device), xb.to(device)
        J_dict, xa, xb, xa2, xb2 = model(xa, xb)
        _, info = loss_fn(J_dict, xa, xb, xa2, xb2, model)

        mixing = _factored_mixing(J_dict, dim_alpha)

        d = dim_alpha + dim_beta
        print(f"\nEvaluation on {num_samples} samples (S^{dim_alpha} × S^{dim_beta}, d={d}):")
        print(f"  Total loss    : {info['total_loss']:.6e}")
        for key in ['nij_00', 'nij_01', 'nij_10', 'nij_11']:
            print(f"  {key}         : {info.get(key, 0.0):.6e}")
        print(f"  Overlap       : {info['overlap_norm']:.6e}")
        print(f"  ACS mixing    : {mixing:.4f}  "
              f"(0 = block-diagonal / product structure, 1 = fully coupled)")

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
        ).to(device)
        model.load_state_dict(checkpoint['model_state_dict'])
        model.eval()

        loss_fn = TotalLoss(manifold_dim=sphere.dim, weights=config['loss_weights']).to(device)

        x = sphere.sample_points_patch1(num_samples).to(device)
        J1, J2, x1, x2 = model(x)
        _, info = loss_fn(J1, J2, x1, x2, model)

        print(f"\nEvaluation on {num_samples} samples (S^{sphere.dim}):")
        print(f"  Total loss : {info['total_loss']:.6e}")
        print(f"  Nij1       : {info['nij1']:.6e}")
        print(f"  Nij2       : {info['nij2']:.6e}")
        print(f"  Overlap    : {info['overlap_norm']:.6e}")

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
