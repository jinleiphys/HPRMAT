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
| 100×100 | 0.004s | 0.000s (27x) | 0.003s (1.3x) | 0.002s (2.0x) | **27x** |
| 400×400 | 0.018s | 0.005s (3.7x) | 0.008s (2.1x) | 0.007s (2.5x) | **3.7x** |
| 1024×1024 | 0.074s | 0.023s (3.2x) | 0.026s (2.8x) | 0.024s (3.1x) | **3.2x** |
| 2000×2000 | 0.226s | 0.077s (2.9x) | 0.070s (3.2x) | 0.061s (3.7x) | **3.7x** |
| 3200×3200 | 0.667s | 0.215s (3.1x) | 0.279s (2.4x) | 0.171s (3.9x) | **3.9x** |
| 4000×4000 | 1.17s | 0.39s (3.0x) | 0.30s (3.9x) | 0.26s (4.4x) | **4.4x** |
| 6400×6400 | 4.45s | 1.19s (3.7x) | 1.22s (3.6x) | 0.83s (5.4x) | **5.4x** |
| 8000×8000 | 7.98s | 2.07s (3.9x) | 2.09s (3.8x) | 1.29s (6.2x) | **6.2x** |
| 10000×10000 | 16.9s | 3.7s (4.6x) | 3.1s (5.5x) | 2.3s (7.4x) | **7.4x** |
| 12800×12800 | 31.4s | 7.4s (4.2x) | 5.3s (5.9x) | 5.2s (6.1x) | **6.1x** |
| 16000×16000 | 60.7s | 14.5s (4.2x) | 9.8s (6.2x) | 8.8s (6.9x) | **6.9x** |
| 25600×25600 | 222.7s | 52.0s (4.3x) | 34.0s (6.6x) | 32.5s (6.9x) | **6.9x** |

#### With GPU (Intel Xeon + RTX 3090)

| Matrix Size | Pierre | Type 1 | Type 2 | Type 3 | Type 4 (GPU) | Best Speedup |
|-------------|--------|--------|--------|--------|--------------|--------------|
| 100×100 | 0.063s | 0.038s (1.7x) | 0.105s (0.6x) | 0.049s (1.3x) | 0.345s (0.2x) | **1.7x** |
| 400×400 | 0.995s | 0.980s (1.0x) | 0.986s (1.0x) | 0.904s (1.1x) | 0.344s (2.9x) | **2.9x** |
| 1024×1024 | 2.32s | 2.27s (1.0x) | 2.45s (0.9x) | 2.37s (1.0x) | 1.12s (2.1x) | **2.1x** |
| 2000×2000 | 3.43s | 3.66s (0.9x) | 3.49s (1.0x) | 3.44s (1.0x) | 1.10s (3.1x) | **3.1x** |
| 4000×4000 | 5.27s | 5.50s (1.0x) | 5.00s (1.1x) | 5.00s (1.1x) | 1.33s (4.0x) | **4.0x** |
| 8000×8000 | 12.5s | 10.7s (1.2x) | 8.19s (1.5x) | 8.38s (1.5x) | 2.01s (6.2x) | **6.2x** |
| 10000×10000 | 18.1s | 10.8s (1.7x) | 8.69s (2.1x) | 9.14s (2.0x) | 2.54s (7.1x) | **7.1x** |
| 12800×12800 | 31.8s | 17.2s (1.8x) | 11.2s (2.8x) | 13.6s (2.3x) | 3.59s (8.8x) | **8.8x** |
| 16000×16000 | 64.3s | 22.3s (2.9x) | 16.4s (3.9x) | 17.7s (3.6x) | 4.72s (13.6x) | **13.6x** |
| 25600×25600 | 216.4s | 104.7s (2.1x) | 41.0s (5.3x) | 41.5s (5.2x) | 12.0s (18.0x) | **18.0x** |

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
