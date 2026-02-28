# Environment Setup

This directory contains the configuration files for setting up the Python environment needed to run the Nijenhuis tensor minimization code.

## Quick Start

### Option 1: Using Conda (Recommended)

1. Create the conda environment:
```bash
conda env create -f environment.yml
```

2. Activate the environment:
```bash
conda activate nijenhuis
```

### Option 2: Using pip with venv

1. Create a virtual environment:
```bash
python -m venv venv
```

2. Activate the virtual environment:
   - On macOS/Linux:
     ```bash
     source venv/bin/activate
     ```
   - On Windows:
     ```bash
     venv\Scripts\activate
     ```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Dependencies

The main dependencies are:

- **PyTorch** (>=2.0.0): Deep learning framework for neural networks and automatic differentiation
- **NumPy** (>=1.24.0): Numerical computing library
- **SciPy** (>=1.10.0): Scientific computing library
- **Matplotlib** (>=3.7.0): Plotting and visualization
- **PyYAML** (>=6.0): YAML configuration file parsing
- **TensorBoard** (>=2.13.0): Training visualization and logging
- **tqdm** (>=4.65.0): Progress bars


## GPU Support

If you have a CUDA-compatible GPU, PyTorch will automatically use it. To verify GPU availability:

```python
import torch
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"CUDA device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A'}")
```

## Updating the Environment

To add new dependencies:

1. For conda: Add them to `environment.yml` and run:
   ```bash
   conda env update -f environment.yml --prune
   ```

2. For pip: Add them to `requirements.txt` and run:
   ```bash
   pip install -r requirements.txt
   ```

## Troubleshooting

### PyTorch Installation Issues

If you encounter issues with PyTorch installation, visit [pytorch.org](https://pytorch.org/get-started/locally/) to get the correct installation command for your system and CUDA version.

### Apple Silicon (M1/M2) Macs

For Apple Silicon Macs, PyTorch has MPS (Metal Performance Shaders) support:

```bash
# Install PyTorch with MPS support
conda install pytorch::pytorch torchvision torchaudio -c pytorch
```

The code will automatically detect and use MPS when available.

## Verifying Installation

After installation, verify everything works:

```bash
python -c "import torch; import numpy; import yaml; print('All imports successful!')"
```

Then test the geometry module:

```bash
cd ..
python -m src.geometry.sphere
```

You should see output indicating successful geometry tests.
