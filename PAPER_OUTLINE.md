# Paper Outline

## Title

**HPRMAT: A high-performance R-matrix solver with GPU acceleration for coupled-channel problems**

## Target Journal

Computer Physics Communications (CPC)

## Authors

[To be added]

## Abstract (Draft)

We present HPRMAT, a high-performance R-matrix solver for coupled-channel scattering problems in nuclear physics. Unlike traditional implementations that rely on matrix inversion, HPRMAT employs direct linear equation solving with optimized BLAS libraries, achieving significant performance improvements. The package provides four solver backends: (1) dense LAPACK using ZGESV, (2) mixed-precision arithmetic with single-precision LU factorization and double-precision iterative refinement, (3) a Woodbury formula approach exploiting the kinetic-coupling matrix structure, and (4) GPU acceleration via NVIDIA cuSOLVER. Benchmark calculations on a 4900×4900 complex matrix (49 channels, 100 Lagrange points) demonstrate that the GPU solver achieves a 3× speedup compared to optimized CPU implementations. The mixed-precision strategy is particularly effective on consumer GPUs where single-precision performance significantly exceeds double-precision. All solvers maintain physics accuracy with relative errors below 10⁻⁵ in cross-section calculations. HPRMAT is designed for modern multi-core CPUs and NVIDIA GPUs, making large-scale coupled-channel calculations accessible on both workstations and high-performance computing clusters.

**Keywords:** R-matrix theory, coupled-channel scattering, GPU acceleration, LAPACK, cuSOLVER, nuclear physics

---

## 1. Introduction

### 1.1 Background
- R-matrix theory in nuclear physics (Lane & Thomas 1958)
- Coupled-channel problems and computational challenges
- Scaling with number of channels: O(N³) for matrix operations

### 1.2 Existing Software
- Descouvemont's R-matrix package (CPC 2016) - Lagrange-Legendre mesh, matrix inversion
- UK PRMAT for atomic physics - MPI parallelization
- FRESCO, CCFULL - different approaches (not R-matrix)

### 1.3 Motivation
- Matrix inversion vs. linear solve (numerical stability, performance)
- GPU acceleration opportunity
- Mixed-precision for consumer GPUs (FP32:FP64 ratio)

### 1.4 Paper Organization
Overview of remaining sections

---

## 2. Theoretical Background

### 2.1 R-Matrix Formulation
- Division of configuration space
- Internal region expansion in Lagrange-Legendre basis
- Boundary conditions and matching

### 2.2 Coupled-Channel Equations
- Channel coupling potentials
- Matrix structure: kinetic + coupling
- Resulting linear system M·x = b

### 2.3 Matrix Structure Analysis
- Block structure: nch × nch blocks of nlag × nlag
- Diagonal blocks: full kinetic matrices
- Off-diagonal blocks: diagonal coupling matrices
- Sparsity analysis and LU fill-in

---

## 3. Numerical Methods

### 3.1 Dense Direct Solve (solver_type=1)
- LAPACK ZGESV implementation
- Reference for accuracy validation

### 3.2 Mixed-Precision Solver (solver_type=2)
- FP32 LU factorization (CGETRF)
- Optional iterative refinement
- Error analysis

### 3.3 Woodbury-Kinetic Method (solver_type=3)
- Matrix decomposition: M = K + V
- Woodbury formula application
- Schur complement approach
- Complexity analysis

### 3.4 GPU Acceleration (solver_type=4)
- cuSOLVER integration
- Host-device memory transfer
- FP64 ↔ FP32 conversion strategy
- Kernel implementation details

---

## 4. Implementation

### 4.1 Code Structure
- Module organization
- Fortran-CUDA interface
- Build system (Makefile, setup.sh)

### 4.2 BLAS Library Integration
- OpenBLAS configuration
- Apple Accelerate comparison
- Thread management

### 4.3 GPU Implementation Details
- Memory allocation strategy
- cuSOLVER handle management
- Error handling

### 4.4 Parallelization Strategy
- Avoiding nested parallelism conflicts
- BLAS threading vs. application-level OpenMP

---

## 5. Performance Results

### 5.1 Test Configuration

**System 1: Apple M3 Ultra (CPU only)**
- CPU: Apple M3 Ultra (32 cores)
- Memory: 96 GB
- BLAS: OpenBLAS (multithreaded)

**System 2: helium1 (CPU + GPU)**
- CPU: Intel Xeon Gold 6248R @ 3.00GHz
- Memory: 630 GB
- GPU: NVIDIA GeForce RTX 3090 (24 GB)
- CUDA: 11.5

### 5.2 CPU-Only Benchmark (Apple M3 Ultra)

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

### 5.3 GPU Benchmark (Intel Xeon + RTX 3090)

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

### 5.4 Solver Comparison Summary

| Solver | Method | Max Speedup | Max Error |
|--------|--------|-------------|-----------|
| Type 1 | Dense LAPACK ZGESV | 5x (CPU) | ~1E-18 |
| Type 2 | Mixed Precision | 7x (CPU) | ~1E-16 |
| Type 3 | Woodbury-Kinetic | 7x (CPU) | ~1E-6 |
| Type 4 | GPU cuSOLVER | **18x** | ~1E-10 |

### 5.5 Key Findings

- **GPU acceleration**: Up to **18x speedup** on RTX 3090 for large matrices
- GPU advantage increases with matrix size (4x at 4000×4000 → 18x at 25600×25600)
- CPU solvers (Type 2, 3) provide 5-7x speedup on Apple M3 Ultra
- All solvers maintain sufficient accuracy for nuclear physics calculations

### 5.6 Accuracy Validation

**Physical Test Cases:**
- Ex1: Alpha-Alpha scattering - S-matrix agrees to 5+ significant digits
- Ex4: 12C+alpha scattering (12 channels) - Amplitudes agree within 0.1%

All solvers validated against Descouvemont's reference code (CPC 200, 2016).

---

## 6. Application Example

### 6.1 Physical Problem
- Deuteron breakup / transfer reaction
- Multi-channel coupling (nch = 49)

### 6.2 Results
- Elastic cross sections
- Timing breakdown
- Practical recommendations

---

## 7. Summary and Future Work

### 7.1 Summary
- Key achievements: first GPU R-matrix solver, 3× speedup
- Solver selection guidelines

### 7.2 Future Directions
- Multi-GPU support
- Batched solves for multiple J values
- Tensor core exploration
- Open-source release plans

---

## Acknowledgments

[Funding, computing resources, etc.]

---

## Appendix A: Installation and Usage

- Prerequisites
- Build instructions
- Input file format
- Example calculations

## Appendix B: Solver Algorithm Pseudocode

Detailed pseudocode for each solver

---

## References

1. A.M. Lane and R.G. Thomas, Rev. Mod. Phys. 30, 257 (1958).
2. D. Baye, Phys. Rep. 565, 1-107 (2015).
3. P. Descouvemont, Comput. Phys. Commun. 200, 199-219 (2016).
4. C.J. Noble et al., Comput. Phys. Commun. 145, 311-340 (2002).
5. E. Anderson et al., LAPACK Users' Guide, 3rd ed., SIAM (1999).
6. NVIDIA cuSOLVER Documentation (2024).
7. [OpenBLAS reference]
8. [Additional physics references for application example]

---

## Figures (Planned)

1. **Fig. 1:** Matrix block structure schematic
2. **Fig. 2:** Sparsity pattern before/after LU factorization
3. **Fig. 3:** Solver comparison bar chart (time)
4. **Fig. 4:** GPU timing breakdown pie chart
5. **Fig. 5:** Scaling with matrix size (log-log plot)
6. **Fig. 6:** Physical results comparison (cross sections)

## Tables (Planned)

1. **Table 1:** Solver type summary
2. **Table 2:** Hardware specifications
3. **Table 3:** Performance benchmark results
4. **Table 4:** Accuracy comparison
5. **Table 5:** Comparison with existing codes

---

## CPC Program Summary

**Title of program:** HPRMAT

**Catalogue identifier:** [To be assigned]

**Program obtainable from:** [GitHub URL / CPC Library]

**Licensing provisions:** [BSD/MIT/GPL]

**Programming language:** Fortran 90/95, CUDA C

**Computer:** Any computer with Fortran compiler; NVIDIA GPU optional

**Operating system:** Linux, macOS

**RAM:** Depends on problem size; ~400 MB for 4900×4900 complex matrix

**Classification:** 17.16 Coupled-channel and distorted-wave methods

**Nature of problem:**
Solving the coupled-channel Schrödinger equation for nuclear scattering using the R-matrix method with Lagrange-Legendre basis functions.

**Solution method:**
Direct solution of the resulting linear system using optimized LAPACK/cuSOLVER routines, with options for mixed-precision and Woodbury formula approaches.

**Restrictions:**
GPU solver requires NVIDIA GPU with CUDA support. Dense storage limits practical matrix sizes to ~10,000.

**Running time:**
Depends on matrix size and hardware. Example: 4900×4900 matrix solves in 0.29s on RTX 3090, 0.89s on 24-core CPU.
