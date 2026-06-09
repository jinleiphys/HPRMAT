!
! gpu_solver_interface.F90 - Fortran interface to GPU solver
!
! Provides a clean Fortran interface to the CUDA solver functions.
! Falls back to CPU solver if GPU is not available.
!
! Usage:
!   use gpu_solver_interface
!   call gpu_solver_init(device_id, ierr)
!   call gpu_solve(A, B, n, nrhs, ierr)
!   call gpu_solver_finalize(ierr)
!

module gpu_solver_interface
  use precision
  use iso_c_binding
  implicit none
  private

  ! Public interface
  public :: gpu_solver_init
  public :: gpu_solver_finalize
  public :: gpu_solve_mixed
  public :: gpu_solve_hybrid
  public :: gpu_solve_double
  public :: gpu_solve_tf32
  public :: gpu_is_available
  public :: gpu_get_info
  public :: gpu_get_device_count
  public :: gpu_multi_init
  public :: gpu_multi_finalize
  public :: gpu_solve_multi
  public :: gpu_solve_multi_mixed
  public :: gpu_solve_auto
  public :: gpu_host_unregister

  ! Module state
  logical, save :: gpu_initialized = .false.
  logical, save :: gpu_available = .false.

  ! C interface declarations
  interface
     ! Initialize GPU solver
     integer(c_int) function gpu_solver_init_c(device_id) bind(C, name="gpu_solver_init_")
       use iso_c_binding
       integer(c_int), intent(in) :: device_id
     end function

     ! Finalize GPU solver
     integer(c_int) function gpu_solver_finalize_c() bind(C, name="gpu_solver_finalize_")
       use iso_c_binding
     end function

     ! Release the cached pinned host-matrix registration (call before freeing the host
     ! matrix buffer on a size change; does not tear down the GPU context)
     integer(c_int) function gpu_host_unregister_c() bind(C, name="gpu_host_unregister_")
       use iso_c_binding
     end function

     ! Mixed precision solver (recommended)
     integer(c_int) function gpu_solve_mixed_c(A, B, n, nrhs, max_refine, tol, info) &
          bind(C, name="gpu_solve_mixed_")
       use iso_c_binding
       complex(c_double_complex), intent(in) :: A(*)
       complex(c_double_complex), intent(inout) :: B(*)
       integer(c_int), intent(in) :: n, nrhs, max_refine
       real(c_double), intent(in) :: tol
       integer(c_int), intent(out) :: info
     end function

     ! Hybrid solver: GPU FP32 factorization + host FP64 residual refinement
     integer(c_int) function gpu_solve_hybrid_c(A, B, n, nrhs, max_refine, tol, info) &
          bind(C, name="gpu_solve_hybrid_")
       use iso_c_binding
       complex(c_double_complex), intent(in) :: A(*)
       complex(c_double_complex), intent(inout) :: B(*)
       integer(c_int), intent(in) :: n, nrhs, max_refine
       real(c_double), intent(in) :: tol
       integer(c_int), intent(out) :: info
     end function

     ! Double precision solver
     integer(c_int) function gpu_solve_double_c(A, B, n, nrhs, info) &
          bind(C, name="gpu_solve_double_")
       use iso_c_binding
       complex(c_double_complex), intent(in) :: A(*)
       complex(c_double_complex), intent(inout) :: B(*)
       integer(c_int), intent(in) :: n, nrhs
       integer(c_int), intent(out) :: info
     end function

     ! TF32 (TensorFloat-32) solver - uses Tensor Core on Ampere+
     integer(c_int) function gpu_solve_tf32_c(A, B, n, nrhs, max_refine, tol, info) &
          bind(C, name="gpu_solve_tf32_")
       use iso_c_binding
       complex(c_double_complex), intent(in) :: A(*)
       complex(c_double_complex), intent(inout) :: B(*)
       integer(c_int), intent(in) :: n, nrhs, max_refine
       real(c_double), intent(in) :: tol
       integer(c_int), intent(out) :: info
     end function

     ! Check GPU availability
     integer(c_int) function gpu_is_available_c() bind(C, name="gpu_is_available_")
       use iso_c_binding
     end function

     ! Get GPU device count
     integer(c_int) function gpu_get_device_count_c() bind(C, name="gpu_get_device_count_")
       use iso_c_binding
     end function

     ! Get GPU info
     subroutine gpu_get_info_c(info_str, len) bind(C, name="gpu_get_info_")
       use iso_c_binding
       character(kind=c_char), intent(out) :: info_str(*)
       integer(c_int), intent(in) :: len
     end subroutine

     ! Multi-GPU initialization
     integer(c_int) function gpu_multi_init_c(ngpu) bind(C, name="gpu_multi_init_")
       use iso_c_binding
       integer(c_int), intent(inout) :: ngpu
     end function

     ! Multi-GPU finalization
     integer(c_int) function gpu_multi_finalize_c() bind(C, name="gpu_multi_finalize_")
       use iso_c_binding
     end function

     ! Multi-GPU solver
     integer(c_int) function gpu_solve_multi_c(A, B, n, nrhs, info) &
          bind(C, name="gpu_solve_multi_")
       use iso_c_binding
       complex(c_double_complex), intent(in) :: A(*)
       complex(c_double_complex), intent(inout) :: B(*)
       integer(c_int), intent(in) :: n, nrhs
       integer(c_int), intent(out) :: info
     end function

     ! Multi-GPU mixed-precision solver (FP32 distributed factor + FP64 host refinement)
     integer(c_int) function gpu_solve_multi_mixed_c(A, B, n, nrhs, max_refine, tol, info) &
          bind(C, name="gpu_solve_multi_mixed_")
       use iso_c_binding
       complex(c_double_complex), intent(in) :: A(*)
       complex(c_double_complex), intent(inout) :: B(*)
       integer(c_int), intent(in) :: n, nrhs, max_refine
       real(c_double), intent(in) :: tol
       integer(c_int), intent(out) :: info
     end function

     ! Auto solver (chooses single/multi GPU automatically)
     integer(c_int) function gpu_solve_auto_c(A, B, n, nrhs, max_refine, tol, info) &
          bind(C, name="gpu_solve_auto_")
       use iso_c_binding
       complex(c_double_complex), intent(in) :: A(*)
       complex(c_double_complex), intent(inout) :: B(*)
       integer(c_int), intent(in) :: n, nrhs, max_refine
       real(c_double), intent(in) :: tol
       integer(c_int), intent(out) :: info
     end function
  end interface

contains

  !---------------------------------------------------------
  ! Initialize GPU solver
  !---------------------------------------------------------
  subroutine gpu_solver_init(device_id, ierr)
    implicit none
    integer, intent(in) :: device_id
    integer, intent(out) :: ierr

    integer(c_int) :: c_device, c_ierr

    ! Check if GPU is available
    if (gpu_is_available_c() == 0) then
       write(*,*) "GPU solver: No GPU available, will use CPU fallback"
       gpu_available = .false.
       ierr = 0
       return
    end if

    gpu_available = .true.
    c_device = device_id
    c_ierr = gpu_solver_init_c(c_device)
    ierr = c_ierr

    if (ierr == 0) then
       gpu_initialized = .true.
    end if
  end subroutine gpu_solver_init

  !---------------------------------------------------------
  ! Finalize GPU solver
  !---------------------------------------------------------
  subroutine gpu_solver_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr

    integer(c_int) :: c_ierr

    if (.not. gpu_initialized) then
       ierr = 0
       return
    end if

    c_ierr = gpu_solver_finalize_c()
    ierr = c_ierr
    gpu_initialized = .false.
  end subroutine gpu_solver_finalize

  ! Release the cached pinned host-matrix registration. Safe to call unconditionally;
  ! it is a no-op if nothing is currently registered.
  subroutine gpu_host_unregister()
    implicit none
    integer(c_int) :: c_ret
    c_ret = gpu_host_unregister_c()
  end subroutine gpu_host_unregister

  !---------------------------------------------------------
  ! Solve using mixed precision (FP32 LU + refinement)
  ! This is the recommended solver for consumer GPUs
  !---------------------------------------------------------
  subroutine gpu_solve_mixed(A, B, n, nrhs, max_refine, tol, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs, max_refine
    real(dp), intent(in) :: tol
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_max_refine, c_info, c_ret
    real(c_double) :: c_tol
    integer :: init_err

    ! Auto-initialize GPU on first call
    if (.not. gpu_initialized) then
       call gpu_solver_init(0, init_err)
    end if

    ! If GPU not available, fall back to CPU
    if (.not. gpu_available) then
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
       return
    end if

    c_n = n
    c_nrhs = nrhs
    c_max_refine = max_refine
    c_tol = tol

    c_ret = gpu_solve_mixed_c(A, B, c_n, c_nrhs, c_max_refine, c_tol, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "GPU solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_mixed

  !---------------------------------------------------------
  ! Hybrid solve: GPU FP32 factorization + host FP64 residual refinement
  !---------------------------------------------------------
  subroutine gpu_solve_hybrid(A, B, n, nrhs, max_refine, tol, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs, max_refine
    real(dp), intent(in) :: tol
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_max_refine, c_info, c_ret
    real(c_double) :: c_tol
    integer :: init_err

    if (.not. gpu_initialized) then
       call gpu_solver_init(0, init_err)
    end if

    if (.not. gpu_available) then
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
       return
    end if

    c_n = n
    c_nrhs = nrhs
    c_max_refine = max_refine
    c_tol = tol

    c_ret = gpu_solve_hybrid_c(A, B, c_n, c_nrhs, c_max_refine, c_tol, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "GPU hybrid solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_hybrid

  !---------------------------------------------------------
  ! Solve using double precision
  !---------------------------------------------------------
  subroutine gpu_solve_double(A, B, n, nrhs, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_info, c_ret

    ! If GPU not available, fall back to CPU
    if (.not. gpu_available) then
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
       return
    end if

    c_n = n
    c_nrhs = nrhs

    c_ret = gpu_solve_double_c(A, B, c_n, c_nrhs, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "GPU solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_double

  !---------------------------------------------------------
  ! Solve using TF32 (TensorFloat-32) precision
  ! Uses Tensor Core acceleration on Ampere+ GPUs (sm_80+)
  ! Falls back to FP32 on older GPUs
  !---------------------------------------------------------
  subroutine gpu_solve_tf32(A, B, n, nrhs, max_refine, tol, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs, max_refine
    real(dp), intent(in) :: tol
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_max_refine, c_info, c_ret
    real(c_double) :: c_tol
    integer :: init_err

    ! Auto-initialize GPU on first call
    if (.not. gpu_initialized) then
       call gpu_solver_init(0, init_err)
    end if

    ! If GPU not available, fall back to CPU
    if (.not. gpu_available) then
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
       return
    end if

    c_n = n
    c_nrhs = nrhs
    c_max_refine = max_refine
    c_tol = tol

    c_ret = gpu_solve_tf32_c(A, B, c_n, c_nrhs, c_max_refine, c_tol, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "TF32 GPU solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_tf32

  !---------------------------------------------------------
  ! Check if GPU is available
  !---------------------------------------------------------
  function gpu_is_available() result(available)
    implicit none
    logical :: available

    available = (gpu_is_available_c() /= 0)
  end function gpu_is_available

  !---------------------------------------------------------
  ! Get GPU info string
  !---------------------------------------------------------
  subroutine gpu_get_info(info_str)
    implicit none
    character(len=*), intent(out) :: info_str

    character(kind=c_char) :: c_str(256)
    integer(c_int) :: c_len
    integer :: i

    c_len = 256
    call gpu_get_info_c(c_str, c_len)

    ! Convert C string to Fortran
    info_str = ''
    do i = 1, 255
       if (c_str(i) == c_null_char) exit
       info_str(i:i) = c_str(i)
    end do
  end subroutine gpu_get_info

  !---------------------------------------------------------
  ! Get number of available GPUs
  !---------------------------------------------------------
  function gpu_get_device_count() result(count)
    implicit none
    integer :: count

    count = gpu_get_device_count_c()
  end function gpu_get_device_count

  !---------------------------------------------------------
  ! Initialize multi-GPU solver
  !---------------------------------------------------------
  subroutine gpu_multi_init(ngpu, ierr)
    implicit none
    integer, intent(inout) :: ngpu
    integer, intent(out) :: ierr

    integer(c_int) :: c_ngpu, c_ret

    c_ngpu = ngpu
    c_ret = gpu_multi_init_c(c_ngpu)
    ngpu = c_ngpu
    ierr = c_ret
  end subroutine gpu_multi_init

  !---------------------------------------------------------
  ! Finalize multi-GPU solver
  !---------------------------------------------------------
  subroutine gpu_multi_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr

    integer(c_int) :: c_ret

    c_ret = gpu_multi_finalize_c()
    ierr = c_ret
  end subroutine gpu_multi_finalize

  !---------------------------------------------------------
  ! Multi-GPU solver
  !---------------------------------------------------------
  subroutine gpu_solve_multi(A, B, n, nrhs, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_info, c_ret

    c_n = n
    c_nrhs = nrhs

    c_ret = gpu_solve_multi_c(A, B, c_n, c_nrhs, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "Multi-GPU solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_multi

  !---------------------------------------------------------
  ! Multi-GPU mixed-precision solver: FP32 distributed LU factorization across the
  ! visible GPUs (cusolverMg) + FP64 host iterative refinement. max_refine controls
  ! the number of refinement steps (0 = bare FP32 solution). Falls back to CPU ZGESV
  ! on any backend error.
  !---------------------------------------------------------
  subroutine gpu_solve_multi_mixed(A, B, n, nrhs, max_refine, tol, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs, max_refine
    real(dp), intent(in) :: tol
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_mr, c_info, c_ret
    real(c_double) :: c_tol

    c_n = n
    c_nrhs = nrhs
    c_mr = max_refine
    c_tol = tol

    c_ret = gpu_solve_multi_mixed_c(A, B, c_n, c_nrhs, c_mr, c_tol, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "Multi-GPU mixed solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_multi_mixed

  !---------------------------------------------------------
  ! Auto solver (automatically chooses single/multi GPU)
  ! This is the recommended solver for most use cases
  !---------------------------------------------------------
  subroutine gpu_solve_auto(A, B, n, nrhs, max_refine, tol, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs, max_refine
    real(dp), intent(in) :: tol
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_max_refine, c_info, c_ret
    real(c_double) :: c_tol
    integer :: init_err

    ! Auto-initialize GPU on first call
    if (.not. gpu_initialized) then
       call gpu_solver_init(0, init_err)
    end if

    ! If GPU not available, fall back to CPU
    if (.not. gpu_available) then
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
       return
    end if

    c_n = n
    c_nrhs = nrhs
    c_max_refine = max_refine
    c_tol = tol

    c_ret = gpu_solve_auto_c(A, B, c_n, c_nrhs, c_max_refine, c_tol, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "GPU solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_auto

  !---------------------------------------------------------
  ! CPU fallback solver (uses LAPACK ZGESV)
  !---------------------------------------------------------
  subroutine cpu_fallback_solve(A, B, n, nrhs, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs
    integer, intent(out) :: ierr

    complex(dp), allocatable :: A_copy(:,:)
    integer, allocatable :: ipiv(:)

    ! Falling back to the CPU: release any cached pinned (page-locked) registration the GPU
    ! backend holds on the host matrix, so a sequence of CPU-fallback solves does not leave
    ! a large page-locked block registered. It is re-pinned on the next successful GPU solve.
    call gpu_host_unregister()

    allocate(A_copy(n, n), ipiv(n))
    A_copy = A

    call zgesv(n, nrhs, A_copy, n, ipiv, B, n, ierr)

    deallocate(A_copy, ipiv)
  end subroutine cpu_fallback_solve

end module gpu_solver_interface
