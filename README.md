# HPRMAT: High-Performance R-Matrix Solver

**A GPU-accelerated R-matrix solver for coupled-channel problems in nuclear physics**

## Overview

HPRMAT is a modern, high-performance implementation of the R-matrix method for solving coupled-channel scattering problems in nuclear physics. It provides significant performance improvements over traditional implementations through:

- **Linear equation solving** instead of matrix inversion
- **GPU acceleration** via NVIDIA cuSOLVER
- **Mixed-precision arithmetic** (FP32 LU factorization + FP64 iterative refinement)
- **TF32 Tensor Core support** for Ampere+ GPUs (RTX 30/40 series)
- **Optimized BLAS libraries** (OpenBLAS with multi-threading)
- **Woodbury formula optimization** exploiting matrix structure
- **Smart setup script** with auto-detection and conda CUDA installation

## Features

| Feature | Description |
|---------|-------------|
| Multiple solver backends | Dense LAPACK, Mixed Precision, Woodbury-Kinetic, GPU cuSOLVER, GPU TF32, Multi-GPU cusolverMg |
| GPU acceleration | ~15x over the optimized CPU direct solver, ~41x over the inversion-based reference (RTX 3090, N=25600, steady state) |
| TF32 Tensor Core | Additional speedup on Ampere+ GPUs (RTX 30/40 series) |
| Parallel BLAS | Full utilization of multi-core CPUs via OpenBLAS |
| Smart setup | Auto-detects GPU, CUDA, OpenBLAS; can install CUDA via conda |
| Lagrange-Legendre basis | Efficient matrix element evaluation with Gauss-Legendre quadrature |
| Open/closed channels | Simultaneous treatment without numerical instability |
| Long-range potentials | Propagation techniques for Coulomb interactions |

## Performance Comparison

### vs. Traditional Matrix Inversion Methods

Traditional R-matrix codes (e.g., Descouvemont 2016) use matrix inversion to solve the linear system. HPRMAT uses LU factorization with forward/backward substitution, which is:

- **More numerically stable** (better condition number handling)
- **Faster** (O(n³) for both, but smaller prefactor for LU solve)
- **Memory efficient** (in-place factorization)

### Benchmark Results (Wall Time)

#### CPU Only (Apple M3 Ultra, 32 cores)

| Matrix Size | Ref. | Type 1 | Type 2 | Type 3 | Best Speedup |
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

#### With GPU (2x Intel Xeon Gold 6248R + RTX 3090, steady state)

Steady-state timings: GPU context and device buffers reused across an energy
scan, pinned host transfer, and the solver reads the caller's matrix directly
(no host copy). CPU columns use fully multi-threaded OpenBLAS
(OPENBLAS_NUM_THREADS=48). These are the numbers published in the CPC paper.

| Matrix Size | Ref. | Type 1 | Type 2 | Type 3 | Type 4 (GPU) | Best Speedup |
|-------------|--------|--------|--------|--------|--------------|--------------|
| 1000×1000 | 1.061s | 1.109s (1.0x) | 1.184s (0.9x) | 1.201s (0.9x) | 0.008s (134x) | **134x** |
| 2000×2000 | 1.508s | 1.516s (1.0x) | 1.681s (0.9x) | 1.770s (0.9x) | 0.021s (70x) | **70x** |
| 4000×4000 | 2.765s | 2.405s (1.1x) | 2.328s (1.2x) | 2.623s (1.1x) | 0.061s (46x) | **46x** |
| 8000×8000 | 9.306s | 5.597s (1.7x) | 4.675s (2.0x) | 5.087s (1.8x) | 0.222s (42x) | **42x** |
| 10000×10000 | 14.409s | 7.347s (2.0x) | 5.739s (2.5x) | 6.515s (2.2x) | 0.365s (40x) | **40x** |
| 16000×16000 | 40.103s | 15.700s (2.6x) | 11.434s (3.5x) | 15.748s (2.5x) | 1.085s (37x) | **37x** |
| 25600×25600 | 144.948s | 52.098s (2.8x) | 41.382s (3.5x) | 48.067s (3.0x) | 3.554s (41x) | **41x** |

Notes: the very large Type 4 speedups at N=1000-2000 are an artifact of
thread-dispatch overhead in the multi-threaded CPU baseline, not a genuine GPU
advantage; the production-regime figure is ~15x over Type 1 (and ~41x over the
reference) at N=25600. When full double precision is required, the host-refined
hybrid mode (`max_refine=2`) takes 9.85s at N=25600 (machine precision, still
5.3x faster than Type 1). The multi-GPU back-end (solver_type=6) extends the
reachable size: N=51200 on one 24 GB card (FP64 matrix kept on the host),
N=64000 on two and N=76800 on four RTX 3090s, in full double precision.

#### CPU Only (Intel i9-12900, 24 threads)

| Matrix Size | Ref. | Type 1 | Type 2 | Type 3 | Best Speedup |
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
| Type 1 | ~1E-18 | Machine precision (matches the reference code) |
| Type 2 | ~1E-15 | Double precision (high accuracy) |
| Type 3 | ~1E-6 | Sufficient for nuclear physics calculations |
| Type 4 | ~1E-10 | GPU FP32 mixed precision (excellent accuracy) |
| Type 5 | ~1E-10 | GPU TF32 Tensor Core (Ampere+ GPUs, slightly faster) |
| Type 6 | ~1E-14 | Multi-GPU (cusolverMg FP32 factorization + FP64 host refinement, default 2 steps); for matrices beyond single-card memory |

### Physical Validation

All solvers validated against Descouvemont's reference code (CPC 200, 2016)
using all five standard test cases from that package:
- Ex1: alpha+208Pb optical model (1 channel) - machine-precision agreement (max deviation ~2E-12 for Type 1)
- Ex2: nucleon-nucleon Reid soft-core (2 channels) - phase shifts agree to better than 1E-4 degrees
- Ex3: 16O+44Ca coupled-channel (4 channels) - cross sections agree within 0.1%
- Ex4: 12C+alpha inelastic (12 channels) - amplitudes agree within 0.1% (1% for Type 3)
- Ex5: non-local Yamaguchi potential (1 channel) - machine-precision phase shifts (Type 1)

## Installation

### Prerequisites

- Fortran compiler (gfortran 9+ recommended)
- LAPACK and BLAS libraries
- OpenBLAS (recommended for CPU performance)
- NVIDIA CUDA Toolkit 11.5+ (for GPU support, 12.0+ for TF32)

### Quick Start (Recommended)

The smart `setup.sh` script automatically detects your environment:

```bash
# Run smart setup (auto-detects GPU, CUDA, OpenBLAS)
./setup.sh

# Build the library
make

# Build examples
make examples
```

The setup script will:
- Detect NVIDIA GPU and determine architecture (sm_86, sm_89, etc.)
- Find CUDA installation (system or conda)
- **Offer to install CUDA 12.9 via conda** if version is too old
- Find OpenBLAS or offer to install it
- Generate optimized `make.inc` configuration

### Installing CUDA via Conda

If you don't have CUDA or need a newer version, the setup script can install it:

```bash
./setup.sh
# Answer 'y' when prompted to install CUDA via conda

# For conda CUDA, activate the environment before building:
source activate_cuda.sh
make
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

Set the module variable `solver_type` (or pass the optional `isolver` argument to `rmatrix`):

```fortran
solver_type = 1   ! Dense LAPACK (ZGESV) - reference implementation
solver_type = 2   ! Mixed Precision (CGETRF + refinement)
solver_type = 3   ! Woodbury-Kinetic - fastest on CPU
solver_type = 4   ! GPU cuSOLVER (FP32) - fastest with NVIDIA GPU
solver_type = 5   ! GPU TF32 - Tensor Core acceleration (Ampere+ GPUs)
solver_type = 6   ! Multi-GPU cusolverMg - for matrices beyond single-card memory
```

The GPU solver (`solver_type = 4`) runs in single precision by default. To recover full
double-precision accuracy on a consumer GPU, pass `max_refine > 0` to enable the
host-refined hybrid mode (FP32 factorization on the GPU, FP64 residual refinement on the
host, where double precision is not throttled):

```fortran
call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, &
                   solver_type=4, max_refine=2)   ! GPU hybrid: FP32 factor + host FP64 refine
```

### Recommended Configuration

| Environment | Recommended Solver |
|-------------|-------------------|
| CPU only (no GPU) | `solver_type = 3` (Woodbury-Kinetic) |
| NVIDIA GPU (any) | `solver_type = 4` (GPU cuSOLVER, FP32) |
| NVIDIA GPU, full double precision | `solver_type = 4, max_refine = 2` (host-refined hybrid) |
| NVIDIA Ampere+ GPU (RTX 30/40) | `solver_type = 5` (GPU TF32) |
| Matrix beyond single-GPU memory | `solver_type = 6` (multi-GPU cusolverMg, FP64 via host refinement) |
| Debugging/validation | `solver_type = 1` (Dense LAPACK) |

### Running the Examples

```bash
make examples          # build all example drivers
cd examples/Ex1
./example1_hp          # alpha+208Pb; prompts for l, nlag, intervals, radius, energies

# or run the full validation sweep over solver types 1-3:
cd examples && ./test_all_solvers.sh
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

#### GPU TF32 (solver_type=5)
Uses TensorFloat-32 format on Ampere+ GPUs (RTX 30/40 series, sm_80+):
1. Enable CUBLAS_TF32_TENSOR_OP_MATH mode
2. LU factorization with Tensor Core acceleration
3. FP64 iterative refinement (5 iterations by default)
4. Falls back to FP32 (Type 4) on older GPUs

TF32 has FP32's dynamic range (8-bit exponent) with reduced mantissa (10-bit), providing ~7-8% speedup over FP32 while maintaining the same accuracy through iterative refinement.

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

### Building and Testing

#### Project Structure

```
HPRMAT/
├── Makefile              # Main library build
├── bindings/
│   ├── Makefile          # Shared library (libhprmat)
│   ├── c/Makefile        # C bindings
│   ├── python/Makefile   # Python bindings (f2py)
│   └── julia/            # Julia module (no build needed)
```

#### Step 1: Build Main Library

```bash
# Configure (auto-detects GPU, CUDA, OpenBLAS)
./setup.sh

# Build static library (libhprmat.a)
make

# Build examples
make examples
```

**Output**: `libhprmat.a` in project root

#### Step 2: Build Shared Library (for bindings)

```bash
cd bindings
make lib
```

**Output**: `lib/libhprmat.dylib` (macOS) or `lib/libhprmat.so` (Linux)

#### Step 3: Build Language Bindings

##### C Bindings

```bash
cd bindings/c
make
```

**Output**: `test_ex1` executable

**Run**:
```bash
# macOS
DYLD_LIBRARY_PATH=../../lib ./test_ex1

# Linux
LD_LIBRARY_PATH=../../lib ./test_ex1
```

##### Python Bindings

```bash
cd bindings/python
make
```

**Output**: `hprmat_fortran.cpython-*.so`

**Test**:
```bash
python -c "from hprmat import RMatrixSolver; print('OK')"
python test_ex1.py
```

##### Julia Bindings

No build needed. Just ensure the shared library is built.

```bash
cd bindings/julia
DYLD_LIBRARY_PATH=../../lib julia test_ex1.jl
```

#### Quick Build (All at Once)

```bash
# From project root
./setup.sh
make
cd bindings && make lib && make python
cd c && make
```

#### Makefile Summary

| Location | Command | Output | Description |
|----------|---------|--------|-------------|
| `/` | `make` | `libhprmat.a` | Static library |
| `/` | `make examples` | `examples/Ex*/` | Example executables |
| `bindings/` | `make lib` | `lib/libhprmat.so` | Shared library |
| `bindings/` | `make python` | Python `.so` | Python f2py module |
| `bindings/c/` | `make` | `test_ex1` | C test executable |
| `bindings/c/` | `make test` | - | Run C test |
| `bindings/python/` | `make` | `hprmat_fortran.so` | Python module |
| `bindings/python/` | `make test` | - | Run Python test |

## Code Structure

```
HPRMAT/
├── src/
│   ├── rmatrix_hp.F90         # Main R-matrix solver interface
│   ├── rmat_solvers.F90       # Six solver implementations (Dense, Mixed, Woodbury, GPU, TF32, Multi-GPU)
│   ├── gpu_solver_interface.F90  # CUDA cuSOLVER wrapper
│   ├── precision.F90          # Precision definitions
│   ├── constants.F90          # Physical constants
│   ├── special_functions.f    # Coulomb/Whittaker functions
│   └── angular_momentum.f     # 3j/6j/9j coefficients
├── examples/
│   ├── Ex0/                   # Benchmarks
│   │   ├── test_solver.f90    # Basic solver test
│   │   ├── benchmark_large.f90  # Performance benchmark
│   │   ├── benchmark_mgpu.f90   # Multi-GPU capacity/scaling benchmark
│   │   ├── benchmark_mgpu_illcond.f90  # Conditioning study (Fig. 5 of the paper)
│   │   └── profile_gpu.f90    # Per-stage GPU profiling (Table 5 of the paper)
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

! Set solver type (1=Dense, 2=Mixed, 3=Woodbury, 4=GPU, 5=TF32, 6=Multi-GPU)
solver_type = 1

! Call R-matrix solver
call rmatrix(nc, lval, qk, eta, rmax, nr, ns, cpot, cu, &
             nmax, nc, nopen, twf, cf, nmax, nc, nc0, nvc, &
             0, cc, solver_type)
```

## Why HPRMAT?

### Comparison with Existing Codes

| Feature | Descouvemont (2016) | UK PRMAT | HPRMAT |
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
- Dense matrix storage (sparse methods not effective due to LU fill-in)
- Multi-GPU support (solver_type = 6, cusolverMg) distributes the FP32
  factorization across visible GPUs and refines in FP64 on the host; its
  benefit is capacity (per-card memory ~1/P), not wall time

## Citation

If you use HPRMAT in your research, please cite:

```bibtex
@article{hprmat2026,
  title={HPRMAT: A high-performance R-matrix solver with GPU acceleration
         for coupled-channel problems in nuclear physics},
  author={Lei, Jin},
  journal={Computer Physics Communications},
  year={2026},
  eprint={2512.11590},
  archivePrefix={arXiv},
  note={in revision}
}
```

## License

MIT License. See the [LICENSE](LICENSE) file for details.

## Authors

Jin Lei

## Acknowledgments

This work was supported by the National Natural Science Foundation of China (Grant Nos. 12475132 and 12535009) and the Fundamental Research Funds for the Central Universities.

---

**HPRMAT** - High-Performance R-Matrix for nuclear coupled-channel calculations
