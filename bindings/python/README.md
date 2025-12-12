# HPRMAT - High-Performance R-Matrix Solver

A GPU-accelerated R-matrix solver for coupled-channel problems in nuclear physics.

## Installation

### From PyPI (Recommended)

```bash
pip install hprmat
```

**Note**: This installs the Python wrapper. You need to build the Fortran extension separately.

### Building from Source

```bash
# Clone the repository
git clone https://github.com/jinlei/hprmat.git
cd hprmat

# Build the library
./setup.sh
make

# Build Python bindings
cd bindings
make python
```

## Quick Start

```python
from hprmat import RMatrixSolver, SOLVER_DENSE
import numpy as np

# Initialize solver
solver = RMatrixSolver(nr=30, ns=1, rmax=10.0)

# Single channel setup
lval = np.array([0], dtype=np.int32)
qk = np.array([0.5], dtype=np.float64)
eta = np.array([0.0], dtype=np.float64)

# Define potential (Gaussian)
cpot = np.zeros((30, 1, 1), dtype=np.complex128, order='F')
for ir, r in enumerate(solver.mesh):
    cpot[ir, 0, 0] = -50.0 * np.exp(-r**2 / 4.0)

# Solve and get S-matrix
S, nopen = solver.solve(lval, qk, eta, cpot)
print(f"S-matrix: {S[0,0]}")
```

## Solver Types

| Type | Method | Best For |
|------|--------|----------|
| `SOLVER_DENSE` (1) | LAPACK ZGESV | Reference, highest precision |
| `SOLVER_MIXED` (2) | Mixed Precision | Large matrices on CPU |
| `SOLVER_WOODBURY` (3) | Woodbury Formula | CPU-only systems |
| `SOLVER_GPU` (4) | NVIDIA cuSOLVER | Systems with GPU |

## Requirements

- Python >= 3.9
- NumPy >= 1.20
- gfortran (for building from source)
- LAPACK/BLAS
- CUDA (optional, for GPU support)

## License

MIT License

## Author

Jin Lei
