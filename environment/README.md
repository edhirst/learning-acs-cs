# Environment Setup

## Setup

```bash
conda create -n nijenhuis python=3.10 -y
conda activate nijenhuis
pip install -r requirements.txt
```

## Dependencies

- **PyTorch** (>=2.0.0)
- **NumPy** (>=1.24.0)
- **SciPy** (>=1.10.0)
- **Matplotlib** (>=3.7.0)
- **PyYAML** (>=6.0)
- **TensorBoard** (>=2.13.0)
- **tqdm** (>=4.65.0)

## GPU Support

The default PyPI `torch` wheel for Linux includes CUDA binaries — no special index URL or CUDA version knowledge is needed. Create the environment on a login node; when the same environment is activated on a GPU node, `torch.cuda.is_available()` returns `True` and the code uses the GPU automatically. On CPU-only machines it falls back to CPU.

For Apple Silicon Macs, MPS acceleration is detected and used automatically.

## Verifying Installation

```bash
python -c "import torch; print('CUDA:', torch.cuda.is_available(), '| MPS:', torch.backends.mps.is_available())"
```
