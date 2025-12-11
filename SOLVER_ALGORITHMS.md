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

**Performance:** ~5% faster than FP64, sufficient accuracy for physics.

### 3.3 Woodbury-Kinetic (solver_type=3)

**Motivation:** Exploit the matrix structure where kinetic energy blocks are separable.

**Matrix decomposition:**
```
M = K_block + V_coupling
```
where K_block is block-diagonal (easily invertible) and V_coupling has low effective rank.

**Woodbury formula:**
```
(A + UCV)⁻¹ = A⁻¹ - A⁻¹U(C⁻¹ + VA⁻¹U)⁻¹VA⁻¹
```

**Algorithm:**
1. Pre-compute K_block⁻¹ for each channel (nch inversions of nlag×nlag)
2. Form Schur complement S = C⁻¹ + V·K_block⁻¹·U
3. Solve S·y = V·K_block⁻¹·b (smaller system)
4. Compute x = K_block⁻¹·b - K_block⁻¹·U·y

**Implementation:**
```fortran
subroutine rmat_cc_woodbury_kinetic(cmat, b, x, n, nch, nlag)
    ! Step 1: Factorize kinetic blocks
    do ich = 1, nch
        call ZGETRF(nlag, nlag, K_inv(ich), nlag, ipiv_k(ich), info)
    end do

    ! Step 2: Build Schur complement (mixed precision)
    ! S_ij = delta_ij - V_ij * K_j^{-1}

    ! Step 3: Solve Schur system
    call CGETRF(nch, nch, S, nch, ipiv_s, info)
    call CGETRS('N', nch, nlag, S, nch, ipiv_s, rhs, nch, info)

    ! Step 4: Back-substitute
    x = K_inv * (b - U * y)
end subroutine
```

**Performance:** ~13% faster than dense ZGESV on CPU.

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
// cusolver_interface.cu
extern "C" void cuda_solve_complex_system_(
    double* A_real, double* A_imag,
    double* b_real, double* b_imag,
    int* n, int* info)
{
    // Allocate device memory
    cuDoubleComplex *d_A, *d_b;
    cuComplex *d_A_f, *d_b_f;  // FP32 versions

    // Copy to device
    cudaMemcpy(d_A, h_A, ...);

    // Convert FP64 -> FP32 on GPU
    convert_z2c_kernel<<<...>>>(d_A, d_A_f, n*n);

    // cuSOLVER LU factorization (FP32)
    cusolverDnCgetrf(handle, n, n, d_A_f, n, workspace, d_ipiv, d_info);

    // Solve (FP32)
    cusolverDnCgetrs(handle, CUBLAS_OP_N, n, 1, d_A_f, n, d_ipiv, d_b_f, n, d_info);

    // Convert FP32 -> FP64 and copy back
    convert_c2z_kernel<<<...>>>(d_b_f, d_b, n);
    cudaMemcpy(h_b, d_b, ...);
}
```

**Fortran Interface:**
```fortran
! gpu_solver_interface.F90
module gpu_solver_interface
    interface
        subroutine cuda_solve_complex_system(A_re, A_im, b_re, b_im, n, info) &
            bind(C, name='cuda_solve_complex_system_')
            use iso_c_binding
            real(c_double) :: A_re(*), A_im(*), b_re(*), b_im(*)
            integer(c_int) :: n, info
        end subroutine
    end interface
end module
```

**Performance:**
| Operation | Time (N=4900) |
|-----------|---------------|
| H2D transfer | 38 ms |
| CGETRF (FP32) | 44 ms |
| CGETRS (FP32) | 1.2 ms |
| D2H transfer | ~1 ms |
| **Total** | **84 ms** |

Compare to CPU ZGESV: 3200 ms → **38x speedup** for the solve kernel.

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

| Solver | Precision | Relative Error | Speed |
|--------|-----------|----------------|-------|
| ZGESV | FP64 | Machine epsilon | 1x |
| Mixed | FP32+refine | ~10⁻¹² | 1.05x |
| Woodbury | FP32 Schur | ~10⁻⁷ | 1.13x |
| GPU | FP32 | ~10⁻⁵ | 3x |

For physics applications, FP32 accuracy (~10⁻⁵ relative error) is typically sufficient.

### 4.3 Iterative Refinement

When higher accuracy is needed with FP32 solvers:

```fortran
! One step of iterative refinement
r = b - matmul(A, x)      ! Residual in FP64
call solve_fp32(A, r, dx) ! Correction in FP32
x = x + dx                ! Update in FP64
```

Typically 1-2 refinement steps recover FP64 accuracy.

## 5. Performance Optimization Tips

### 5.1 BLAS Library Selection

**Recommendation: Use OpenBLAS on all platforms (including Apple Silicon)**

Our benchmarks show OpenBLAS outperforms Apple Accelerate for complex LAPACK routines:

| Platform | BLAS Library | ZGESV Time (4900×4900) | Relative |
|----------|--------------|------------------------|----------|
| Apple M3 Max | **OpenBLAS** | **27.1 s** | **1.0x** |
| Apple M3 Max | Apple Accelerate | 50.0 s | 0.54x |

Apple Accelerate saturates at ~4 threads for complex matrix operations, while OpenBLAS scales effectively to all cores.

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

For repeated solves, keep data on GPU:

```fortran
! Initialize once
call cuda_init_solver(n)

! Repeated solves (data stays on GPU)
do iter = 1, niter
    call cuda_solve_inplace(A, b, x)
end do

! Cleanup
call cuda_finalize_solver()
```

## 6. Future Directions

### 6.1 Multi-GPU Support

Distribute different J (total angular momentum) values across multiple GPUs:

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
