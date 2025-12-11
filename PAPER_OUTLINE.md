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

*[Section to be completed with systematic benchmark data]*

### 5.1 Test Configuration
- Hardware specifications (to be documented)
- Test cases: varying nch (10, 25, 49), nlag (50, 100, 150)
- Metrics: wall-clock time, FLOPS, memory usage

### 5.2 Comparison with Descouvemont's Code
**[KEY COMPARISON - Required for publication]**
- Same physical problem solved by both codes
- Wall-clock time comparison
- Matrix inversion (Descouvemont) vs. linear solve (HPRMAT)
- Scaling with matrix size

### 5.3 HPRMAT Solver Comparison
- solver_type 1 (Dense LAPACK) - baseline
- solver_type 2 (Mixed Precision)
- solver_type 3 (Woodbury-Kinetic)
- solver_type 4 (GPU cuSOLVER)

### 5.4 GPU vs CPU Performance
- RTX 3090 vs OpenBLAS multi-thread
- Time breakdown: transfer, factorization, solve
- Scaling with matrix size

### 5.5 Accuracy Validation
- Cross-section comparison across all solvers
- Comparison with Descouvemont's results
- Relative error analysis

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
