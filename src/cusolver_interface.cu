/*
 * cusolver_interface.cu - Optimized GPU solver for PSTARS
 *
 * Uses FP32 (single precision) for LU factorization on GPU
 * with optional FP64 iterative refinement for accuracy.
 *
 * Key optimizations:
 * - FP64→FP32 conversion on GPU using CUDA kernels
 * - Persistent memory allocation (reused across calls)
 * - Minimized CPU-GPU data transfers
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cuComplex.h>

// Error checking macros
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(err)); \
        return -1; \
    } \
} while(0)

#define CUSOLVER_CHECK(call) do { \
    cusolverStatus_t status = call; \
    if (status != CUSOLVER_STATUS_SUCCESS) { \
        fprintf(stderr, "cuSOLVER error at %s:%d: %d\n", __FILE__, __LINE__, status); \
        return -1; \
    } \
} while(0)

#define CUBLAS_CHECK(call) do { \
    cublasStatus_t status = call; \
    if (status != CUBLAS_STATUS_SUCCESS) { \
        fprintf(stderr, "cuBLAS error at %s:%d: %d\n", __FILE__, __LINE__, status); \
        return -1; \
    } \
} while(0)

// Global handles
static cusolverDnHandle_t cusolver_handle = NULL;
static cublasHandle_t cublas_handle = NULL;
static int gpu_initialized = 0;
static int current_device = 0;

// Persistent GPU memory
static cuComplex *d_A_sp = NULL;       // Single precision matrix
static cuComplex *d_B_sp = NULL;       // Single precision RHS
static cuDoubleComplex *d_A_dp = NULL; // Double precision matrix (for refinement)
static cuDoubleComplex *d_B_dp = NULL; // Double precision RHS
static cuDoubleComplex *d_R = NULL;    // Residual
static cuDoubleComplex *d_X = NULL;    // Solution
static int *d_ipiv = NULL;
static int *d_info = NULL;
static void *d_work = NULL;
static size_t work_size = 0;
static int alloc_n = 0;
static int alloc_nrhs = 0;

// ============================================
// CUDA Kernels for FP64 <-> FP32 conversion
// ============================================

// Convert complex double to complex float (on GPU)
__global__ void convert_z2c_kernel(const cuDoubleComplex* __restrict__ src,
                                    cuComplex* __restrict__ dst,
                                    int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = make_cuFloatComplex((float)src[idx].x, (float)src[idx].y);
    }
}

// Convert complex float to complex double (on GPU)
__global__ void convert_c2z_kernel(const cuComplex* __restrict__ src,
                                    cuDoubleComplex* __restrict__ dst,
                                    int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = make_cuDoubleComplex((double)src[idx].x, (double)src[idx].y);
    }
}

// Compute residual R = B - A*X on GPU
__global__ void compute_residual_kernel(const cuDoubleComplex* __restrict__ A,
                                         const cuDoubleComplex* __restrict__ X,
                                         const cuDoubleComplex* __restrict__ B,
                                         cuDoubleComplex* __restrict__ R,
                                         int n, int nrhs) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    int col = blockIdx.y;

    if (row < n && col < nrhs) {
        cuDoubleComplex sum = make_cuDoubleComplex(0.0, 0.0);
        for (int k = 0; k < n; k++) {
            cuDoubleComplex a = A[row + k * n];
            cuDoubleComplex x = X[k + col * n];
            sum.x += a.x * x.x - a.y * x.y;
            sum.y += a.x * x.y + a.y * x.x;
        }
        R[row + col * n].x = B[row + col * n].x - sum.x;
        R[row + col * n].y = B[row + col * n].y - sum.y;
    }
}

// Add correction: X += dX
__global__ void add_correction_kernel(cuDoubleComplex* __restrict__ X,
                                       const cuComplex* __restrict__ dX_sp,
                                       int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        X[idx].x += (double)dX_sp[idx].x;
        X[idx].y += (double)dX_sp[idx].y;
    }
}

// ============================================
// Initialization and cleanup
// ============================================

extern "C" {

int gpu_solver_init_(int *device_id) {
    if (gpu_initialized) {
        return 0;
    }

    current_device = *device_id;
    CUDA_CHECK(cudaSetDevice(current_device));

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, current_device));
    printf("GPU Solver initialized on: %s (%.1f GB)\n",
           prop.name, prop.totalGlobalMem / 1e9);

    CUSOLVER_CHECK(cusolverDnCreate(&cusolver_handle));
    CUBLAS_CHECK(cublasCreate(&cublas_handle));

    gpu_initialized = 1;
    return 0;
}

int gpu_solver_finalize_() {
    if (!gpu_initialized) {
        return 0;
    }

    if (d_A_sp) cudaFree(d_A_sp);
    if (d_B_sp) cudaFree(d_B_sp);
    if (d_A_dp) cudaFree(d_A_dp);
    if (d_B_dp) cudaFree(d_B_dp);
    if (d_R) cudaFree(d_R);
    if (d_X) cudaFree(d_X);
    if (d_ipiv) cudaFree(d_ipiv);
    if (d_info) cudaFree(d_info);
    if (d_work) cudaFree(d_work);

    d_A_sp = d_B_sp = NULL;
    d_A_dp = d_B_dp = d_R = d_X = NULL;
    d_ipiv = d_info = NULL;
    d_work = NULL;
    alloc_n = alloc_nrhs = 0;

    if (cusolver_handle) cusolverDnDestroy(cusolver_handle);
    if (cublas_handle) cublasDestroy(cublas_handle);

    cusolver_handle = NULL;
    cublas_handle = NULL;
    gpu_initialized = 0;

    return 0;
}

static int ensure_gpu_memory(int n, int nrhs) {
    if (n <= alloc_n && nrhs <= alloc_nrhs) {
        return 0;
    }

    // Free old allocations
    if (d_A_sp) cudaFree(d_A_sp);
    if (d_B_sp) cudaFree(d_B_sp);
    if (d_A_dp) cudaFree(d_A_dp);
    if (d_B_dp) cudaFree(d_B_dp);
    if (d_R) cudaFree(d_R);
    if (d_X) cudaFree(d_X);
    if (d_ipiv) cudaFree(d_ipiv);
    if (d_info) cudaFree(d_info);
    if (d_work) cudaFree(d_work);

    // Allocate new memory
    CUDA_CHECK(cudaMalloc(&d_A_sp, n * n * sizeof(cuComplex)));
    CUDA_CHECK(cudaMalloc(&d_B_sp, n * nrhs * sizeof(cuComplex)));
    CUDA_CHECK(cudaMalloc(&d_A_dp, n * n * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_B_dp, n * nrhs * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_R, n * nrhs * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_X, n * nrhs * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_ipiv, n * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_info, sizeof(int)));

    // Get workspace size for CGETRF
    int lwork;
    CUSOLVER_CHECK(cusolverDnCgetrf_bufferSize(cusolver_handle, n, n, d_A_sp, n, &lwork));

    work_size = lwork * sizeof(cuComplex);
    CUDA_CHECK(cudaMalloc(&d_work, work_size));

    alloc_n = n;
    alloc_nrhs = nrhs;

    return 0;
}

// ============================================
// Optimized Mixed Precision Solver
// ============================================

int gpu_solve_mixed_(
    double *h_A,
    double *h_B,
    int *n_ptr,
    int *nrhs_ptr,
    int *max_refine_ptr,
    double *tol_ptr,
    int *info_ptr
) {
    int n = *n_ptr;
    int nrhs = *nrhs_ptr;
    int max_refine = *max_refine_ptr;
    double tol = *tol_ptr;
    int h_info;

    if (!gpu_initialized) {
        int device = 0;
        if (gpu_solver_init_(&device) != 0) {
            *info_ptr = -1;
            return -1;
        }
    }

    if (ensure_gpu_memory(n, nrhs) != 0) {
        *info_ptr = -2;
        return -1;
    }

    cuDoubleComplex *h_A_dp = (cuDoubleComplex*)h_A;
    cuDoubleComplex *h_B_dp = (cuDoubleComplex*)h_B;

    // Copy FP64 data to GPU
    CUDA_CHECK(cudaMemcpy(d_A_dp, h_A_dp, n * n * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_B_dp, h_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));

    // Convert FP64 -> FP32 on GPU (fast!)
    int blockSize = 256;
    int numBlocks_A = (n * n + blockSize - 1) / blockSize;
    int numBlocks_B = (n * nrhs + blockSize - 1) / blockSize;

    convert_z2c_kernel<<<numBlocks_A, blockSize>>>(d_A_dp, d_A_sp, n * n);
    convert_z2c_kernel<<<numBlocks_B, blockSize>>>(d_B_dp, d_B_sp, n * nrhs);
    CUDA_CHECK(cudaGetLastError());

    // LU factorization in FP32
    CUSOLVER_CHECK(cusolverDnCgetrf(cusolver_handle, n, n, d_A_sp, n,
                                     (cuComplex*)d_work, d_ipiv, d_info));
    CUDA_CHECK(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        fprintf(stderr, "CGETRF failed with info = %d\n", h_info);
        *info_ptr = h_info;
        return -1;
    }

    // Solve in FP32
    CUSOLVER_CHECK(cusolverDnCgetrs(cusolver_handle, CUBLAS_OP_N, n, nrhs,
                                     d_A_sp, n, d_ipiv, d_B_sp, n, d_info));

    // Convert solution FP32 -> FP64 on GPU
    convert_c2z_kernel<<<numBlocks_B, blockSize>>>(d_B_sp, d_X, n * nrhs);
    CUDA_CHECK(cudaGetLastError());

    // Iterative refinement on GPU (if requested)
    if (max_refine > 0) {
        cuDoubleComplex alpha = make_cuDoubleComplex(-1.0, 0.0);
        cuDoubleComplex beta = make_cuDoubleComplex(1.0, 0.0);

        for (int iter = 0; iter < max_refine; iter++) {
            // Compute residual R = B - A*X on GPU using cuBLAS
            CUDA_CHECK(cudaMemcpy(d_R, d_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                                  cudaMemcpyDeviceToDevice));
            CUBLAS_CHECK(cublasZgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                      n, nrhs, n,
                                      &alpha, d_A_dp, n, d_X, n,
                                      &beta, d_R, n));

            // Check convergence (compute max norm on GPU)
            int max_idx;
            CUBLAS_CHECK(cublasIzamax(cublas_handle, n * nrhs, d_R, 1, &max_idx));
            cuDoubleComplex max_val;
            CUDA_CHECK(cudaMemcpy(&max_val, d_R + max_idx - 1, sizeof(cuDoubleComplex),
                                  cudaMemcpyDeviceToHost));
            double max_res = sqrt(max_val.x * max_val.x + max_val.y * max_val.y);

            if (max_res < tol) {
                break;
            }

            // Convert residual to FP32 and solve for correction
            convert_z2c_kernel<<<numBlocks_B, blockSize>>>(d_R, d_B_sp, n * nrhs);
            CUSOLVER_CHECK(cusolverDnCgetrs(cusolver_handle, CUBLAS_OP_N, n, nrhs,
                                             d_A_sp, n, d_ipiv, d_B_sp, n, d_info));

            // Add correction: X += dX (on GPU)
            add_correction_kernel<<<numBlocks_B, blockSize>>>(d_X, d_B_sp, n * nrhs);
            CUDA_CHECK(cudaGetLastError());
        }
    }

    // Copy solution back to host
    CUDA_CHECK(cudaMemcpy(h_B_dp, d_X, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyDeviceToHost));

    *info_ptr = 0;
    return 0;
}

// Double precision solver (for comparison)
int gpu_solve_double_(
    double *h_A,
    double *h_B,
    int *n_ptr,
    int *nrhs_ptr,
    int *info_ptr
) {
    int n = *n_ptr;
    int nrhs = *nrhs_ptr;
    int h_info;

    if (!gpu_initialized) {
        int device = 0;
        if (gpu_solver_init_(&device) != 0) {
            *info_ptr = -1;
            return -1;
        }
    }

    if (ensure_gpu_memory(n, nrhs) != 0) {
        *info_ptr = -2;
        return -1;
    }

    cuDoubleComplex *h_A_dp = (cuDoubleComplex*)h_A;
    cuDoubleComplex *h_B_dp = (cuDoubleComplex*)h_B;

    CUDA_CHECK(cudaMemcpy(d_A_dp, h_A_dp, n * n * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_B_dp, h_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyHostToDevice));

    // Get workspace for ZGETRF
    int lwork;
    CUSOLVER_CHECK(cusolverDnZgetrf_bufferSize(cusolver_handle, n, n, d_A_dp, n, &lwork));

    if (lwork * sizeof(cuDoubleComplex) > work_size) {
        if (d_work) cudaFree(d_work);
        work_size = lwork * sizeof(cuDoubleComplex);
        CUDA_CHECK(cudaMalloc(&d_work, work_size));
    }

    CUSOLVER_CHECK(cusolverDnZgetrf(cusolver_handle, n, n, d_A_dp, n,
                                     (cuDoubleComplex*)d_work, d_ipiv, d_info));
    CUDA_CHECK(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        *info_ptr = h_info;
        return -1;
    }

    CUSOLVER_CHECK(cusolverDnZgetrs(cusolver_handle, CUBLAS_OP_N, n, nrhs,
                                     d_A_dp, n, d_ipiv, d_B_dp, n, d_info));

    CUDA_CHECK(cudaMemcpy(h_B_dp, d_B_dp, n * nrhs * sizeof(cuDoubleComplex),
                          cudaMemcpyDeviceToHost));

    *info_ptr = 0;
    return 0;
}

int gpu_is_available_() {
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        return 0;
    }
    return 1;
}

void gpu_get_info_(char *info_str, int *len) {
    if (!gpu_is_available_()) {
        snprintf(info_str, *len, "No GPU available");
        return;
    }

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    snprintf(info_str, *len, "%s (%.1f GB)", prop.name, prop.totalGlobalMem / 1e9);
}

} // extern "C"
