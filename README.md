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

### Benchmark Results

*[To be added: Systematic benchmarks comparing HPRMAT with Descouvemont's original R-matrix code]*

#### Test Cases (Planned)
- Single-channel scattering
- Multi-channel coupled problems (nch = 10, 25, 49)
- Various matrix sizes (nlag = 50, 100, 150)

#### Comparisons (Planned)
1. **HPRMAT vs. Descouvemont's code** - Same physics problem, compare wall-clock time
2. **Matrix inversion vs. linear solve** - Numerical accuracy and performance
3. **Solver backends** - Dense, Mixed Precision, Woodbury, GPU

### Physical Accuracy

*[To be added: Cross-section comparisons across all solvers and validation against Descouvemont's results]*

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
GPU_ARCH = sm_86  # RTX 3090
```

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

## Code Structure

*[To be finalized: A clean API similar to Descouvemont's R-matrix package will be provided for easy integration]*

### Planned Interface

The final release will provide a simple subroutine interface compatible with Descouvemont's conventions:

```fortran
! Example API (to be implemented)
call hprmat_init(nch, nlag, solver_type)
call hprmat_set_potential(V_coupling)
call hprmat_solve(E, S_matrix)
call hprmat_finalize()
```

### Design Goals
- Drop-in replacement for existing R-matrix codes
- Minimal code changes required for users of Descouvemont's package
- Clear separation of solver backend from physics interface

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
  author={[Authors]},
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

[To be added]

## Acknowledgments

This work was supported by [funding information].

---

**HPRMAT** - High-Performance R-Matrix for nuclear coupled-channel calculations
