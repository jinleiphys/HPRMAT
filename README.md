# HPRMAT: High-Performance R-Matrix Solver

**A GPU-accelerated R-matrix solver for coupled-channel problems in nuclear physics**

## Overview

HPRMAT is a modern, high-performance implementation of the R-matrix method for solving coupled-channel scattering problems in nuclear physics. It provides significant performance improvements over traditional implementations through:

- **Linear equation solving** instead of matrix inversion
- **GPU acceleration** via NVIDIA cuSOLVER
- **Mixed-precision arithmetic** (FP32 LU factorization + FP64 iterative refinement)
- **Optimized BLAS libraries** (OpenBLAS with multi-threading)
- **Woodbury formula optimization** exploiting matrix structure

## Features

| Feature | Description |
|---------|-------------|
| Multiple solver backends | Dense LAPACK, Mixed Precision, Woodbury-Kinetic, GPU cuSOLVER |
| GPU acceleration | 3x speedup using NVIDIA cuSOLVER with FP32 |
| Parallel BLAS | Full utilization of multi-core CPUs via OpenBLAS |
| Lagrange-Legendre basis | Efficient matrix element evaluation with Gauss-Legendre quadrature |
| Open/closed channels | Simultaneous treatment without numerical instability |
| Long-range potentials | Propagation techniques for Coulomb interactions |

## Performance Comparison

### vs. Traditional Matrix Inversion Methods

Traditional R-matrix codes (e.g., Descouvemont 2015) use matrix inversion to solve the linear system. HPRMAT uses LU factorization with forward/backward substitution, which is:

- **More numerically stable** (better condition number handling)
- **Faster** (O(n³) for both, but smaller prefactor for LU solve)
- **Memory efficient** (in-place factorization)

### Benchmark Results (Wall Time)

#### CPU Only (Apple M3 Ultra, 32 cores)

| Matrix Size | Pierre | Type 1 | Type 2 | Type 3 | Best Speedup |
|-------------|--------|--------|--------|--------|--------------|
| 100×100 | 0.002s | 0.000s (17x) | 0.002s (1.0x) | 0.002s (1.3x) | **17x** |
| 400×400 | 0.018s | 0.005s (3.6x) | 0.009s (2.1x) | 0.007s (2.4x) | **3.6x** |
| 1024×1024 | 0.091s | 0.024s (3.7x) | 0.036s (2.5x) | 0.025s (3.6x) | **3.7x** |
| 2000×2000 | 0.489s | 0.077s (6.3x) | 0.069s (7.1x) | 0.060s (8.2x) | **8.2x** |
| 3200×3200 | 1.357s | 0.208s (6.5x) | 0.303s (4.5x) | 0.173s (7.8x) | **7.8x** |
| 4000×4000 | 1.198s | 0.364s (3.3x) | 0.338s (3.5x) | 0.280s (4.3x) | **4.3x** |
| 6400×6400 | 4.334s | 1.185s (3.7x) | 1.258s (3.4x) | 0.818s (5.3x) | **5.3x** |
| 8000×8000 | 8.162s | 2.066s (4.0x) | 1.965s (4.2x) | 1.319s (6.2x) | **6.2x** |
| 10000×10000 | 15.6s | 3.68s (4.2x) | 3.05s (5.1x) | 2.31s (6.8x) | **6.8x** |
| 12800×12800 | 31.7s | 7.41s (4.3x) | 5.50s (5.8x) | 5.06s (6.3x) | **6.3x** |
| 16000×16000 | 60.5s | 14.1s (4.3x) | 10.3s (5.9x) | 8.82s (6.9x) | **6.9x** |
| 25600×25600 | 231.5s | 52.3s (4.4x) | 34.0s (6.8x) | 32.0s (7.2x) | **7.2x** |

#### With GPU (Intel Xeon + RTX 3090)

| Matrix Size | Pierre | Type 1 | Type 2 | Type 3 | Type 4 (GPU) | Best Speedup |
|-------------|--------|--------|--------|--------|--------------|--------------|
| 100×100 | 0.055s | 0.036s (1.5x) | 0.113s (0.5x) | 0.050s (1.1x) | 0.364s (0.2x) | **1.5x** |
| 400×400 | 1.102s | 1.110s (1.0x) | 1.042s (1.1x) | 0.941s (1.2x) | 0.349s (3.2x) | **3.2x** |
| 1024×1024 | 2.403s | 2.363s (1.0x) | 2.586s (0.9x) | 2.647s (0.9x) | 1.050s (2.3x) | **2.3x** |
| 2000×2000 | 3.466s | 3.555s (1.0x) | 3.528s (1.0x) | 3.575s (1.0x) | 1.107s (3.1x) | **3.1x** |
| 4000×4000 | 5.196s | 4.787s (1.1x) | 5.206s (1.0x) | 4.926s (1.1x) | 1.311s (4.0x) | **4.0x** |
| 8000×8000 | 12.4s | 9.69s (1.3x) | 7.65s (1.6x) | 7.39s (1.7x) | 2.07s (6.0x) | **6.0x** |
| 10000×10000 | 17.4s | 10.6s (1.6x) | 9.13s (1.9x) | 8.98s (1.9x) | 2.68s (6.5x) | **6.5x** |
| 12800×12800 | 34.4s | 15.4s (2.2x) | 11.7s (2.9x) | 13.1s (2.6x) | 3.68s (9.3x) | **9.3x** |
| 16000×16000 | 62.9s | 21.7s (2.9x) | 15.2s (4.1x) | 17.3s (3.6x) | 4.99s (12.6x) | **12.6x** |
| 25600×25600 | 210.0s | 105.2s (2.0x) | 40.9s (5.1x) | 39.6s (5.3x) | 11.8s (17.8x) | **17.8x** |

#### CPU Only (Intel i9-12900, 24 threads)

| Matrix Size | Pierre | Type 1 | Type 2 | Type 3 | Best Speedup |
|-------------|--------|--------|--------|--------|--------------|
| 100×100 | 0.058s | 0.000s (331x) | 0.053s (1.1x) | 0.000s (301x) | **331x** |
| 400×400 | 0.077s | 0.003s (31x) | 0.025s (3.1x) | 0.002s (38x) | **38x** |
| 1024×1024 | 0.146s | 0.227s (0.6x) | 0.020s (7.4x) | 0.021s (6.9x) | **7.4x** |
| 2000×2000 | 0.294s | 0.393s (0.8x) | 0.123s (2.4x) | 0.114s (2.6x) | **2.6x** |
| 3200×3200 | 1.115s | 0.557s (2.0x) | 0.414s (2.7x) | 0.336s (3.3x) | **3.3x** |
| 4000×4000 | 2.207s | 0.840s (2.6x) | 0.539s (4.1x) | 0.555s (4.0x) | **4.1x** |
| 6400×6400 | 8.44s | 3.00s (2.8x) | 2.79s (3.0x) | 2.36s (3.6x) | **3.6x** |
| 8000×8000 | 21.9s | 7.04s (3.1x) | 6.51s (3.4x) | 4.76s (4.6x) | **4.6x** |
| 10000×10000 | 48.1s | 18.5s (2.6x) | 12.1s (4.0x) | 11.3s (4.3x) | **4.3x** |
| 12800×12800 | 99.5s | 35.8s (2.8x) | 22.9s (4.4x) | 23.4s (4.3x) | **4.4x** |
| 16000×16000 | 163.0s | 54.7s (3.0x) | 36.9s (4.4x) | 32.9s (5.0x) | **5.0x** |
| 25600×25600 | 679.3s | 217.6s (3.1x) | 125.1s (5.4x) | 120.1s (5.7x) | **5.7x** |

### Solver Accuracy

| Solver | Max Error | Description |
|--------|-----------|-------------|
| Type 1 | ~1E-18 | Machine precision (identical to Pierre's original) |
| Type 2 | ~1E-16 | Double precision (high accuracy) |
| Type 3 | ~1E-6 | Sufficient for nuclear physics calculations |
| Type 4 | ~1E-10 | GPU mixed precision (excellent accuracy) |

### Physical Validation

All solvers validated against Descouvemont's reference code (CPC 200, 2016):
- Ex1: Alpha-Alpha scattering (1 channel) - Phase shifts match to 5+ digits
- Ex4: 12C+alpha scattering (12 channels) - Amplitudes agree within 0.1%

## Installation

### Prerequisites

- Fortran compiler (gfortran 9+ recommended)
- LAPACK and BLAS libraries
- OpenBLAS (recommended for CPU performance)
- NVIDIA CUDA Toolkit 11.5+ (for GPU support)

### Building

```bash
# Configure build (auto-detects GPU and OpenBLAS)
./setup.sh

# Build the solver
cd cookie && make clean && make
```

### Manual Configuration

Edit `make.inc` to specify:

```makefile
# CPU-only build (Mac/Linux)
OPENBLAS_DIR = /opt/homebrew/opt/openblas
LIBSTD1 = -L$(OPENBLAS_DIR)/lib -lopenblas -fopenmp

# GPU-enabled build (Linux with NVIDIA GPU)
GPU_ENABLED = true
CUDA_PATH = /usr/local/cuda
NVCC = $(CUDA_PATH)/bin/nvcc
GPU_ARCH = sm_86  # Set according to your GPU (see table below)
```

#### GPU Architecture Values

| GPU | Architecture | GPU_ARCH |
|-----|--------------|----------|
| RTX 2080/2070/2060 | Turing | sm_75 |
| RTX 3090/3080/3070/3060 | Ampere | sm_86 |
| RTX 4090/4080/4070/4060 | Ada Lovelace | sm_89 |
| RTX 5090/5080/5070/5060 | Blackwell | sm_120 |
| V100 | Volta | sm_70 |
| A100 | Ampere | sm_80 |
| H100 | Hopper | sm_90 |

Note: The `setup.sh` script auto-detects your GPU architecture. Manual configuration is only needed if auto-detection fails.

## Usage

### Solver Selection

Set `solver_type` in your input file:

```fortran
solver_type = 1   ! Dense LAPACK (ZGESV) - reference implementation
solver_type = 2   ! Mixed Precision (CGETRF + refinement)
solver_type = 3   ! Woodbury-Kinetic - fastest on CPU
solver_type = 4   ! GPU cuSOLVER - fastest with NVIDIA GPU
```

### Recommended Configuration

| Environment | Recommended Solver |
|-------------|-------------------|
| CPU only (no GPU) | `solver_type = 3` (Woodbury-Kinetic) |
| NVIDIA GPU available | `solver_type = 4` (GPU cuSOLVER) |
| Debugging/validation | `solver_type = 1` (Dense LAPACK) |

### Example Input

```
test_calculation.in
---
! R-matrix coupled-channel calculation
solver_type = 4        ! Use GPU solver
nlag = 100             ! Lagrange points per channel
...
```

### Running

```bash
cd cookie
echo "test/test_multibin_dense.in" | ./cookie
```

## Theory

### R-Matrix Method

The R-matrix method divides configuration space into an internal region (r < a) and external region (r > a). In the internal region, the wave function is expanded in a complete basis set (Lagrange-Legendre functions based on Gauss-Legendre quadrature), leading to a linear system:

```
(H - E) Ψ = 0  →  M · X = B
```

where M is the Hamiltonian matrix with boundary conditions.

### Matrix Structure

For nch channels and nlag Lagrange points, the matrix M has dimension (nch × nlag) × (nch × nlag) with block structure:

```
      ch1           ch2           ch3
    [  K_full   ] [  Diag    ] [  Diag    ]  ch1
    [  Diag     ] [  K_full  ] [  Diag    ]  ch2
    [  Diag     ] [  Diag    ] [  K_full  ]  ch3
```

- **Diagonal blocks**: Full kinetic energy matrix (nlag × nlag)
- **Off-diagonal blocks**: Diagonal coupling potential matrices

### Solver Algorithms

#### Dense LAPACK (solver_type=1)
Standard LU factorization using ZGESV. Reference implementation.

#### Mixed Precision (solver_type=2)
1. Convert matrix to single precision (complex64)
2. LU factorization using CGETRF (FP32)
3. Solve using CGETRS (FP32)
4. Optional: iterative refinement in FP64

#### Woodbury-Kinetic (solver_type=3)
Exploits the matrix structure using the Woodbury formula:
```
(K + UV)⁻¹ = K⁻¹ - K⁻¹U(I + VK⁻¹U)⁻¹VK⁻¹
```
where K is block-diagonal (kinetic energy) and UV represents coupling.

#### GPU cuSOLVER (solver_type=4)
1. Transfer matrix to GPU memory
2. Convert FP64 → FP32 on GPU
3. LU factorization using cuSOLVER CGETRF
4. Solve using CGETRS
5. Transfer solution back to CPU

## Language Bindings

HPRMAT provides bindings for multiple programming languages:

| Language | Method | Directory |
|----------|--------|-----------|
| **C/C++** | ISO_C_BINDING | `bindings/c/` |
| **Python** | f2py (direct Fortran) | `bindings/python/` |
| **Julia** | ccall (direct Fortran) | `bindings/julia/` |

### Python Example

```python
from hprmat import RMatrixSolver, SOLVER_DENSE
import numpy as np

# Initialize solver
solver = RMatrixSolver(nr=30, ns=1, rmax=10.0)

# Single channel setup
lval = np.array([0], dtype=np.int32)
qk = np.array([0.5], dtype=np.float64)
eta = np.array([0.0], dtype=np.float64)

# Define potential
cpot = np.zeros((30, 1, 1), dtype=np.complex128, order='F')
for ir, r in enumerate(solver.mesh):
    cpot[ir, 0, 0] = -50.0 * np.exp(-r**2 / 4.0)

# Solve and get S-matrix
S, nopen = solver.solve(lval, qk, eta, cpot)
print(f"S-matrix: {S[0,0]}")
```

### Julia Example

```julia
include("HPRMAT.jl")
using .HPRMAT

# Initialize solver
solver = RMatrixSolver(30, 1, 10.0)

# Single channel setup
lval = Int32[0]
qk = [0.5]
eta = [0.0]

# Define potential
cpot = zeros(ComplexF64, 30, 1, 1)
for (ir, r) in enumerate(solver.mesh)
    cpot[ir, 1, 1] = -50.0 * exp(-r^2 / 4.0)
end

# Solve and get S-matrix
S, nopen = solve(solver, lval, qk, eta, cpot)
println("S-matrix: ", S[1,1])
```

### Building Bindings

```bash
cd bindings
make lib       # Build shared library
make python    # Build Python module
```

See `bindings/README.md` for detailed documentation.

## Code Structure

```
HPRMAT/
├── src/
│   ├── rmatrix_hp.F90         # Main R-matrix solver interface
│   ├── rmat_solvers.F90       # Four solver implementations (Dense, Mixed, Woodbury, GPU)
│   ├── gpu_solver_interface.F90  # CUDA cuSOLVER wrapper
│   ├── precision.F90          # Precision definitions
│   ├── constants.F90          # Physical constants
│   ├── special_functions.f    # Coulomb/Whittaker functions
│   └── angular_momentum.f     # 3j/6j/9j coefficients
├── examples/
│   ├── Ex0/                   # Large matrix benchmark
│   │   ├── test_solver.f90    # Basic solver test
│   │   └── benchmark_large.f90  # Performance benchmark
│   ├── Ex1/example1_hp.f90    # Alpha-Alpha scattering (1 channel)
│   ├── Ex2/example2_hp.f90    # Reid NN potential (2 channels)
│   ├── Ex3/example3_hp.f90    # 16O+44Ca scattering (4 channels)
│   ├── Ex4/example4_hp.f90    # 12C+alpha scattering (12 channels)
│   └── Ex5/example5_hp.f90    # Yamaguchi non-local potential
├── bindings/
│   ├── c/                     # C/C++ interface
│   ├── python/                # Python wrapper
│   └── julia/                 # Julia module
├── rmat_pierre/               # Pierre Descouvemont's original code (reference)
└── lib/                       # Compiled libraries
```

### Fortran API

```fortran
use rmat_hp_mod

! Set solver type (1=Dense, 2=Mixed, 3=Woodbury, 4=GPU)
solver_type = 1

! Call R-matrix solver
call rmatrix(nc, lval, qk, eta, rmax, nr, ns, cpot, cu, &
             nmax, nc, nopen, twf, cf, nmax, nc, nc0, nvc, &
             0, cc, solver_type)
```

## Why HPRMAT?

### Comparison with Existing Codes

| Feature | Descouvemont (2015) | UK PRMAT | HPRMAT |
|---------|--------------------|---------:|--------|
| Target | Nuclear physics | Atomic physics | Nuclear physics |
| Parallelization | None | MPI + OpenMP | OpenBLAS + GPU |
| GPU support | No | No | **Yes (cuSOLVER)** |
| Linear solve | Matrix inversion | ScaLAPACK | LAPACK/cuSOLVER |
| Mixed precision | No | No | **Yes** |

### Key Innovations

1. **First GPU-accelerated R-matrix solver for nuclear physics**
2. **Mixed-precision strategy** optimized for RTX consumer GPUs (FP32:FP64 = 64:1)
3. **Woodbury formula exploitation** of kinetic-coupling matrix structure
4. **Modern BLAS integration** (OpenBLAS, cuBLAS)

## Limitations

- GPU solver requires NVIDIA GPU with CUDA support
- Single-GPU implementation (multi-GPU planned)
- Dense matrix storage (sparse methods not effective due to LU fill-in)

## Citation

If you use HPRMAT in your research, please cite:

```bibtex
@article{hprmat2025,
  title={HPRMAT: A high-performance R-matrix solver with GPU acceleration
         for coupled-channel problems in nuclear physics},
  author={Lei, Jin},
  journal={Computer Physics Communications},
  year={2025},
  note={In preparation}
}
```

## References

1. P. Descouvemont, "An R-matrix package for coupled-channel problems in nuclear physics,"
   Comput. Phys. Commun. 200, 199-219 (2016). [arXiv:1510.03540](https://arxiv.org/abs/1510.03540)

2. D. Baye, "The Lagrange-mesh method," Phys. Rep. 565, 1-107 (2015).

3. A.M. Lane and R.G. Thomas, "R-Matrix Theory of Nuclear Reactions,"
   Rev. Mod. Phys. 30, 257 (1958).

4. C.J. Noble et al., "A parallel R-matrix program PRMAT for electron-atom
   and electron-ion scattering calculations," Comput. Phys. Commun. 145, 311-340 (2002).

## License

[To be determined]

## Authors

Jin Lei (jinleiphys@gmail.com)

## Acknowledgments

This work was supported by [funding information].

---

**HPRMAT** - High-Performance R-Matrix for nuclear coupled-channel calculations
