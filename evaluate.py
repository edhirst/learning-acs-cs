"""
Utility script to evaluate a trained model checkpoint.
"""

import torch
import argparse

from src.geometry import EvenSphere
from src.models import TwoPatchAlmostComplexNet
from src.loss import TotalLoss


def evaluate_model(checkpoint_path: str, num_samples: int = 1000):
    """
    Evaluate a trained model and print loss metrics.

    Args:
        checkpoint_path: Path to model checkpoint.
        num_samples: Number of collocation points for evaluation.
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    config = checkpoint['config']

    print(f"Checkpoint epoch : {checkpoint['epoch']}")
    print(f"Best training loss: {checkpoint['best_loss']:.6e}")

    sphere = EvenSphere(dimension=config['geometry']['dimension'], device=device)

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
