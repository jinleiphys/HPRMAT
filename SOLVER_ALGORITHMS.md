# HPRMAT Solver Algorithms

This document provides detailed technical descriptions of the solver algorithms implemented in HPRMAT.

## 1. Problem Formulation

### 1.1 Coupled-Channel Scattering

In nuclear scattering with coupling between channels, we solve the coupled Schrödinger equations:

```
[-ℏ²/(2μ) d²/dr² + l(l+1)ℏ²/(2μr²) + V_αα(r) - E] ψ_α(r) + Σ_α' V_αα'(r) ψ_α'(r) = 0
```

where α labels the channels (including internal quantum numbers and angular momentum).

### 1.2 R-Matrix Formulation

The R-matrix method divides space at radius r = a:

**Internal region (r < a):** Expand in Lagrange-Legendre basis functions (using Gauss-Legendre quadrature points):
```
ψ_α(r) = Σ_i c_αi f_i(r)
```

**External region (r > a):** Match to Coulomb/asymptotic wave functions.

### 1.3 Linear System

The expansion coefficients satisfy:

```
M · c = b
```

where:
- M is the Hamiltonian matrix of dimension N = nch × nlag
- c is the coefficient vector
- b is the boundary condition vector

## 2. Matrix Structure Analysis

### 2.1 Block Structure

The matrix M has nch × nch blocks, each of size nlag × nlag:

```
M = | K₁ + D₁₁   D₁₂       D₁₃      ...  |
    | D₂₁        K₂ + D₂₂  D₂₃      ...  |
    | D₃₁        D₃₂       K₃ + D₃₃ ...  |
    | ...        ...       ...       ... |
```

where:
- K_α = Kinetic energy matrix (FULL nlag × nlag matrix)
- D_αα' = Coupling potential matrix (DIAGONAL in Lagrange basis)

### 2.2 Sparsity Analysis

For typical calculations (nlag=100, nch=49):
- Total elements: 4900² = 24,010,000
- Diagonal blocks: 49 × 100² = 490,000 (full)
- Off-diagonal blocks: 49 × 48 × 100 = 235,200 (diagonal only)
- **Initial sparsity: ~3%**

**Critical observation:** Despite initial sparsity, LU factorization causes fill-in, making sparse solvers ineffective.

### 2.3 Why Sparse Solvers Fail

We tested several sparse approaches:
- UMFPACK: Slower than dense for N > 2000 due to fill-in overhead
- ILU-GMRES: 50x slower than direct solve
- Block Gauss-Seidel: Does not converge (strong coupling)

**Conclusion:** Dense direct solvers are optimal for this problem structure.

## 3. Solver Implementations

### 3.1 Dense LAPACK (solver_type=1)

**Algorithm:**
```fortran
call ZGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
```

**Complexity:** O(N³) for factorization, O(N²) for solve

**Advantages:**
- Reference implementation
- Maximum numerical stability
- Well-tested LAPACK routines

**Implementation:**
```fortran
subroutine rmat_cc_origin(cmat, b, x, n, nch, nlag)
    ! Direct solve using LAPACK ZGESV
    complex(8), intent(inout) :: cmat(n,n)  ! Matrix (overwritten)
    complex(8), intent(inout) :: b(n)       ! RHS (overwritten with solution)

    call ZGESV(n, 1, cmat, n, ipiv, b, n, info)
    x = b
end subroutine
```

### 3.2 Mixed Precision (solver_type=2)

**Motivation:** Single precision (FP32) operations are 2x faster than double precision (FP64) on most hardware, and much faster on GPUs.

**Algorithm:**
1. Convert A (complex128) → A_s (complex64)
2. LU factorization: A_s = P·L·U using CGETRF
3. Solve: L·U·x_s = P·b using CGETRS
4. Convert x_s (complex64) → x (complex128)
5. Optional: Iterative refinement in FP64

**Implementation:**
```fortran
subroutine rmat_cc_mixed_precision(cmat, b, x, n, nch, nlag)
    complex(8) :: cmat(n,n), b(n), x(n)
    complex(4), allocatable :: cmat_s(:,:), b_s(:), x_s(:)

    ! Convert to single precision
    cmat_s = cmplx(cmat, kind=4)
    b_s = cmplx(b, kind=4)

    ! Solve in single precision
    call CGETRF(n, n, cmat_s, n, ipiv, info)
    call CGETRS('N', n, 1, cmat_s, n, ipiv, b_s, n, info)

    ! Convert back
    x = cmplx(b_s, kind=8)

    ! Optional: iterative refinement
    ! r = b - A*x
    ! Solve A*dx = r
    ! x = x + dx
end subroutine
```

**Performance:** 5-7x faster than the inversion-based reference for large matrices (N > 10000), with ~1E-15 accuracy.

### 3.3 Woodbury-Kinetic (solver_type=3)

**Motivation:** Exploit the block structure: the kinetic block K is real and
shared across all channels, so K⁻¹ is computed only once, in fast real
arithmetic.

**Matrix decomposition:**
```
M = K_diag + V_coupling = K_diag (I + K_diag⁻¹ V_coupling)
```
where K_diag is block-diagonal (kinetic + diagonal potentials) and V_coupling
holds the off-diagonal channel coupling.

**Important:** V_coupling is NOT low-rank, so the literal Woodbury identity
`(A + UCV)⁻¹ = A⁻¹ - A⁻¹U(C⁻¹ + VA⁻¹U)⁻¹VA⁻¹` by itself gives no cost
reduction here; it is quoted in the paper only as motivating background. The
implemented algorithm solves the Schur-type system `(I + K⁻¹V) c = K⁻¹ b`,
which is still N×N: the dominant cost remains an O(N³) single-precision
factorization, and the gain is a constant factor (real arithmetic for K⁻¹,
single precision for the Schur solve), not a change in asymptotic complexity.

**Algorithm:**
1. Extract and factorize the shared real kinetic block K once (DGETRF, ~(2/3)·nlag³ real ops)
2. Extract the nch×nch coupling matrices V(r_i) at each Lagrange point
3. Build the N×N Schur matrix S, S[(i-1)nch+a, (j-1)nch+b] = δ_ij δ_ab + K⁻¹(i,j) V_ab(r_j)
4. Transform the right-hand side: b̃ = K⁻¹ b
5. Solve S·x̃ = b̃ in single precision (CGETRF/CGETRS), transform back

**Implementation sketch:**
```fortran
! Step 1: factorize the shared real kinetic block once
call DGETRF(nlag, nlag, K, nlag, ipiv_k, info)

! Steps 2-3: build the N x N Schur matrix S = I + K^{-1} V

! Step 5: solve the Schur system in single precision
call CGETRF(ntot, ntot, S, ntot, ipiv_s, info)
call CGETRS('N', ntot, nch, S, ntot, ipiv_s, rhs, ntot, info)
```

**Performance:** Up to 7.2x faster than the inversion-based reference for large
matrices (Apple M3 Ultra, N=25600), with ~1E-6 accuracy (sufficient for nuclear
physics; not recommended near narrow resonances).

### 3.4 GPU cuSOLVER (solver_type=4)

**Motivation:** Modern GPUs offer massive parallelism for dense linear algebra. RTX 3090 has 35 TFLOPS FP32 vs 0.56 TFLOPS FP64.

**Strategy:** Use FP32 on GPU (where FP32:FP64 = 64:1 performance ratio).

**Algorithm:**
1. Allocate GPU memory
2. Copy matrix A and vector b to GPU (Host → Device)
3. Convert FP64 → FP32 on GPU (CUDA kernel)
4. LU factorization using cuSOLVER CGETRF
5. Solve using cuSOLVER CGETRS
6. Convert FP32 → FP64 on GPU
7. Copy solution x back to CPU (Device → Host)

**CUDA Interface:**
```c
// cusolver_interface.cu (actual entry point; matrix and RHS are passed as
// double-precision complex arrays, all conversion happens on the device)
extern "C" int gpu_solve_mixed_(
    double *h_A,        /* complex(dp) matrix, column-major */
    double *h_B,        /* complex(dp) RHS, overwritten with solution */
    int *n_ptr, int *nrhs_ptr, int *max_refine_ptr,
    double *tol_ptr, int *info_ptr)
{
    // Persistent device buffers (allocated on first call, reused after)
    cuDoubleComplex *d_A, *d_B;
    cuComplex *d_A_f, *d_B_f;  // FP32 versions

    // Copy to device (pinned host registration is cached)
    cudaMemcpy(d_A, h_A, ...);

    // Convert FP64 -> FP32 on GPU
    convert_z2c_kernel<<<...>>>(d_A, d_A_f, n*n);

    // cuSOLVER LU factorization (FP32)
    cusolverDnCgetrf(handle, n, n, d_A_f, n, workspace, d_ipiv, d_info);

    // Solve (FP32)
    cusolverDnCgetrs(handle, CUBLAS_OP_N, n, nrhs, d_A_f, n, d_ipiv, d_B_f, n, d_info);

    // Convert FP32 -> FP64 and copy back
    convert_c2z_kernel<<<...>>>(d_B_f, d_B, n*nrhs);
    cudaMemcpy(h_B, d_B, ...);
}
```

**Fortran Interface:**
```fortran
! gpu_solver_interface.F90
interface
   integer(c_int) function gpu_solve_mixed_c(A, B, n, nrhs, max_refine, tol, info) &
        bind(C, name="gpu_solve_mixed_")
     use iso_c_binding
     complex(c_double_complex), intent(in)    :: A(*)
     complex(c_double_complex), intent(inout) :: B(*)
     integer(c_int), intent(in) :: n, nrhs, max_refine
     real(c_double), intent(in) :: tol
     integer(c_int), intent(out) :: info
   end function
end interface
! Analogous bindings: gpu_solve_hybrid_ (host-refined), gpu_solve_multi_mixed_ (multi-GPU)
```

**Performance (per-stage, N=25600, RTX 3090, steady state, pinned transfer):**

| Stage | Time |
|-------|------|
| H2D transfer (pinned, 10.5 GB) | 0.86 s |
| FP64→FP32 conversion | 0.02 s |
| FP32 LU factorization (cusolverDnCgetrf) | 2.49 s |
| FP32 triangular solves (cusolverDnCgetrs) | 0.08 s |
| FP32→FP64 conversion + D2H | 0.02 s |
| Host RHS assembly + R-matrix extraction | 0.07 s |
| **Total (steady state)** | **3.54 s** |

Compare 52.1 s for the multi-threaded CPU direct solver (Type 1) and 144.9 s
for the inversion-based reference: ~15x and ~41x respectively. The one-time GPU
context creation and host page-locking (~1.1 s) are paid once per energy scan.

## 4. Numerical Considerations

### 4.1 Condition Number

The coupled-channel matrix can be ill-conditioned for:
- Near-threshold energies
- Strong coupling potentials
- Many closed channels

**Monitoring:**
```fortran
! Estimate condition number
call ZGECON('1', n, A, n, anorm, rcond, work, rwork, info)
if (rcond < 1e-12) warning("Ill-conditioned matrix")
```

### 4.2 Accuracy vs Speed Trade-off

**Benchmark Results (Wall Time, Large Matrices)**

| Matrix Size | Ref. (inversion) | Type 1 | Type 2 | Type 3 | Best Speedup |
|-------------|-------------------|--------|--------|--------|--------------|
| 4000×4000 | 1.198s | 0.364s (3.3x) | 0.338s (3.5x) | 0.280s (4.3x) | **4.3x** |
| 10000×10000 | 15.6s | 3.68s (4.2x) | 3.05s (5.1x) | 2.31s (6.8x) | **6.8x** |
| 25600×25600 | 231.5s | 52.3s (4.4x) | 34.0s (6.8x) | 32.0s (7.2x) | **7.2x** |

(Apple M3 Ultra, 32 cores; these are the published CPC numbers.)

**Accuracy Comparison**

| Solver | Max Error | Description |
|--------|-----------|-------------|
| Type 1 (ZGESV) | ~1E-18 | Machine precision |
| Type 2 (Mixed) | ~1E-15 | Double precision |
| Type 3 (Woodbury) | ~1E-6 | Sufficient for nuclear physics |

For nuclear physics applications, Type 3 accuracy (~10⁻⁶ relative error) is sufficient while providing the best performance.

### 4.3 Iterative Refinement

When higher accuracy is needed with FP32 solvers:

```fortran
! One step of iterative refinement
r = b - matmul(A, x)      ! Residual in FP64
call solve_fp32(A, r, dx) ! Correction in FP32
x = x + dx                ! Update in FP64
```

Typically 1-2 refinement steps recover FP64 accuracy on well-conditioned
systems. Convergence requires the single-precision factor to be a contraction
(κ·u_FP32 < 1, i.e. κ below roughly 1E7); beyond that the refinement cannot
converge and the solver falls back to a full double-precision solve (see the
conditioning study, Fig. 5 of the paper).

## 5. Performance Optimization Tips

### 5.1 BLAS Library Selection

**Recommendation: Use OpenBLAS on all platforms (including Apple Silicon)**

Our benchmarks on **Apple M3 Ultra (32 cores, 96GB RAM)** show OpenBLAS provides excellent multithreaded performance:

| Matrix Size | Ref. (inversion) | Type 1 (ZGESV) | Type 3 (Woodbury) |
|-------------|-------------------|----------------|-------------------|
| 10000×10000 | 15.6s | 3.68s | 2.31s |
| 25600×25600 | 231.5s | 52.3s | 32.0s |

OpenBLAS scales effectively to all cores for complex matrix operations.

```bash
# OpenBLAS (recommended)
export OPENBLAS_NUM_THREADS=24

# Avoid nested parallelism when using BLAS threading
export OMP_NUM_THREADS=1
```

**Installation:**
```bash
# macOS
brew install openblas

# Linux
apt install libopenblas-dev  # or from source
```

### 5.2 Memory Layout

Fortran uses column-major storage. Ensure matrices are accessed column-wise:

```fortran
! Good: column-wise access
do j = 1, n
    do i = 1, n
        A(i,j) = ...
    end do
end do

! Bad: row-wise access (cache unfriendly)
do i = 1, n
    do j = 1, n
        A(i,j) = ...
    end do
end do
```

### 5.3 GPU Memory Management

Device buffers, cuSOLVER handles, and the pinned host-matrix registration are
persistent: they are created on the first call and reused across an energy scan.
No explicit init call is needed:

```fortran
use gpu_solver_interface

! Repeated solves (buffers and handles reused automatically)
do iter = 1, niter
    call gpu_solve_mixed(A, X, n, nrhs, max_refine, tol, ierr)
end do

! Cleanup at end of run (also releases the pinned registration)
call gpu_solver_finalize()
```

## 6. Multi-GPU Support and Future Directions

### 6.1 Multi-GPU Support (implemented, solver_type = 6)

For matrices beyond single-card memory, `gpu_solve_multi_mixed_` distributes the
FP32 LU factorization across the visible GPUs with cusolverMg (1D block-cyclic
column layout) and recovers full double precision through FP64 iterative
refinement on the host (default 2 steps). Per-card device memory scales as
roughly 1/P; the benefit is capacity, not wall time. See the manuscript section
"Scalability and multi-GPU support" and `examples/Ex0/benchmark_mgpu.f90`.

A complementary future direction is distributing different J (total angular
momentum) values across GPUs, one independent system per card:

```
GPU 0: J = 0, 4, 8, ...
GPU 1: J = 1, 5, 9, ...
GPU 2: J = 2, 6, 10, ...
GPU 3: J = 3, 7, 11, ...
```

### 6.2 Batched Solves

Solve multiple systems simultaneously using cuSOLVER batched routines:

```c
cusolverDnCgetrfBatched(handle, n, d_Aarray, n, d_ipivArray, d_infoArray, batchSize);
```

### 6.3 Tensor Core Acceleration

For supported GPUs (Volta+), explore FP16 with tensor cores:
- NVIDIA A100: 312 TFLOPS FP16 tensor
- Requires careful numerical analysis

---

## References

1. E. Anderson et al., "LAPACK Users' Guide," 3rd ed., SIAM, 1999.
2. NVIDIA, "cuSOLVER Library Documentation," 2024.
3. M. Frigo and S.G. Johnson, "The Design and Implementation of FFTW3," Proc. IEEE 93, 216-231 (2005).
4. J. Dongarra et al., "The International Exascale Software Project Roadmap," Int. J. High Perf. Comp. App. 25, 3-60 (2011).
