# HPRMAT Language Bindings

This directory contains language bindings for HPRMAT (High-Performance R-Matrix Solver).

## Supported Languages

| Language | Method | Directory |
|----------|--------|-----------|
| **C/C++** | ISO_C_BINDING | `c/` |
| **Python** | f2py (direct Fortran) | `python/` |
| **Julia** | ccall (direct Fortran) | `julia/` |

## Quick Start

### Build Shared Library

```bash
cd bindings
make lib
```

This creates `lib/libhprmat.so` (Linux) or `lib/libhprmat.dylib` (macOS).

### Python

```bash
cd bindings/python
make
python example_alpha_alpha.py
```

```python
from hprmat import RMatrixSolver, SOLVER_DENSE
import numpy as np

# Initialize
solver = RMatrixSolver(nr=30, ns=1, rmax=10.0)

# Single channel
lval = np.array([0], dtype=np.int32)
qk = np.array([0.5], dtype=np.float64)
eta = np.array([0.0], dtype=np.float64)

# Potential
cpot = np.zeros((30, 1, 1), dtype=np.complex128, order='F')
for ir, r in enumerate(solver.mesh):
    cpot[ir, 0, 0] = -50.0 * np.exp(-r**2 / 4.0)

# Solve
S, nopen = solver.solve(lval, qk, eta, cpot)
print(f"S-matrix: {S[0,0]}")
```

### Julia

```bash
cd bindings
make lib  # Build shared library first

cd julia
julia example_alpha_alpha.jl
```

```julia
include("HPRMAT.jl")
using .HPRMAT

# Initialize
solver = RMatrixSolver(30, 1, 10.0)

# Single channel
lval = Int32[0]
qk = [0.5]
eta = [0.0]

# Potential
cpot = zeros(ComplexF64, 30, 1, 1)
for (ir, r) in enumerate(solver.mesh)
    cpot[ir, 1, 1] = -50.0 * exp(-r^2 / 4.0)
end

# Solve
S, nopen = solve(solver, lval, qk, eta, cpot)
println("S-matrix: ", S[1,1])
```

### C/C++

```bash
cd bindings
make lib  # Build shared library first

cd c
make
./example_alpha_alpha
```

```c
#include "hprmat.h"

// Initialize
double zrma[30];
hprmat_init(30, 1, 10.0, zrma);

// Set solver
hprmat_set_solver(HPRMAT_SOLVER_DENSE);

// Solve
int nopen;
double _Complex cu[1];
hprmat_solve(nch, lval, qk, eta, rmax, nr, ns, cpot, nr*ns, cu, &nopen, 0);
```

## Solver Types

| Solver | Constant | Description |
|--------|----------|-------------|
| 1 | `SOLVER_DENSE` | Dense LAPACK ZGESV (reference, most accurate) |
| 2 | `SOLVER_MIXED` | Mixed precision (faster for large matrices) |
| 3 | `SOLVER_WOODBURY` | CPU optimized (5-7x speedup) |
| 4 | `SOLVER_GPU` | GPU cuSOLVER (requires CUDA, up to 18x speedup) |

## API Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `hprmat_init(nr, ns, rmax, zrma)` | Initialize R-matrix calculation |
| `hprmat_solve(...)` | Solve scattering problem |
| `hprmat_solve_full(...)` | Solve with wave function output |
| `hprmat_wavefunction(...)` | Interpolate wave function |
| `hprmat_set_solver(type)` | Set solver type |

### Parameters

- `nr`: Number of Lagrange functions per interval (typically 20-40)
- `ns`: Number of intervals (typically 1)
- `rmax`: R-matrix channel radius (fm)
- `nch`: Number of channels
- `lval`: Angular momentum for each channel
- `qk`: Wave numbers (positive=open, negative=closed)
- `eta`: Sommerfeld parameters (for Coulomb)
- `cpot`: Coupling potential matrix

## File Structure

```
bindings/
├── Makefile              # Main build system
├── README.md             # This file
├── c/
│   ├── hprmat.h          # C header
│   ├── hprmat_c.F90      # ISO_C_BINDING interface
│   ├── example_alpha_alpha.c
│   └── Makefile
├── python/
│   ├── hprmat.py         # Python wrapper class
│   ├── hprmat_f2py.F90   # f2py interface
│   ├── example_alpha_alpha.py
│   └── Makefile
└── julia/
    ├── HPRMAT.jl         # Julia module
    └── example_alpha_alpha.jl
```

## Requirements

- **Fortran**: gfortran with OpenMP support
- **LAPACK/BLAS**: System libraries or OpenBLAS
- **Python**: NumPy, f2py (part of NumPy)
- **Julia**: Julia 1.6+
- **C**: Any C99 compiler

## Notes

- Python and Julia wrappers call Fortran directly (no C layer)
- C interface is provided for C/C++, Rust, Go, etc.
- All arrays use Fortran (column-major) memory layout
- Complex numbers use double precision (complex*16 / ComplexF64)
