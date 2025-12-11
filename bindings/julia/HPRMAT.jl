"""
HPRMAT.jl - High-Performance R-Matrix Solver for Julia

Julia interface for the HPRMAT Fortran library using ccall.

# Example
```julia
using HPRMAT

# Initialize
nr, ns, rmax = 20, 1, 10.0
solver = RMatrixSolver(nr, ns, rmax)

# Single channel problem
nch = 1
lval = Int32[0]
qk = [0.5]
eta = [0.0]

# Build potential
cpot = zeros(ComplexF64, nr*ns, nch, nch)
for (ir, r) in enumerate(solver.mesh)
    cpot[ir, 1, 1] = -50.0 * exp(-r^2 / 4.0)
end

# Solve
S, nopen = solve(solver, lval, qk, eta, cpot)
println("S-matrix: ", S[1,1])
```

Author: Jin Lei
Date: December 2025
"""
module HPRMAT

using LinearAlgebra

export RMatrixSolver, solve, solve_with_wavefunction, solve_nonlocal
export interpolate_wavefunction, set_solver!
export SOLVER_DENSE, SOLVER_MIXED, SOLVER_WOODBURY, SOLVER_GPU

# Solver type constants
const SOLVER_DENSE = 1      # Dense LAPACK ZGESV (reference)
const SOLVER_MIXED = 2      # Mixed precision
const SOLVER_WOODBURY = 3   # Woodbury-Kinetic (CPU optimized)
const SOLVER_GPU = 4        # GPU cuSOLVER

const SOLVER_NAMES = Dict(
    SOLVER_DENSE => "Dense LAPACK ZGESV (reference)",
    SOLVER_MIXED => "Mixed Precision (single + double refinement)",
    SOLVER_WOODBURY => "Woodbury-Kinetic (CPU optimized)",
    SOLVER_GPU => "GPU cuSOLVER (NVIDIA GPU)",
)

# Library path - will be set when module loads
const libhprmat = Ref{String}("")

function __init__()
    # Search for the library
    search_paths = [
        joinpath(@__DIR__, "..", "..", "lib"),
        joinpath(@__DIR__, "..", "lib"),
        "/usr/local/lib",
        "/usr/lib",
        get(ENV, "HPRMAT_LIB_PATH", ""),
    ]

    lib_names = Sys.isapple() ? ["libhprmat.dylib", "libhprmat.so"] :
                Sys.iswindows() ? ["hprmat.dll", "libhprmat.dll"] :
                ["libhprmat.so"]

    for path in search_paths
        isempty(path) && continue
        for name in lib_names
            lib_path = joinpath(path, name)
            if isfile(lib_path)
                libhprmat[] = lib_path
                return
            end
        end
    end

    @warn "HPRMAT library not found. Build it first with 'make lib'"
end

"""
Check if library is loaded.
"""
function check_lib()
    if isempty(libhprmat[])
        error("HPRMAT library not found. Set HPRMAT_LIB_PATH or build with 'make lib'")
    end
end

"""
    RMatrixSolver

High-performance R-matrix solver for nuclear scattering.

# Fields
- `nr::Int`: Number of Lagrange functions per interval
- `ns::Int`: Number of intervals
- `rmax::Float64`: R-matrix channel radius (fm)
- `solver::Int`: Solver type (1-4)
- `mesh::Vector{Float64}`: Lagrange mesh abscissas
"""
mutable struct RMatrixSolver
    nr::Int
    ns::Int
    rmax::Float64
    solver::Int
    mesh::Vector{Float64}

    function RMatrixSolver(nr::Integer, ns::Integer, rmax::Real;
                           solver::Integer=SOLVER_DENSE)
        check_lib()

        nr = Int32(nr)
        ns = Int32(ns)
        rmax = Float64(rmax)
        mesh = zeros(Float64, nr * ns)

        # Call Fortran: subroutine rmat_ini(nr, ns, rmax, zrma)
        ccall((:rmat_ini_, libhprmat[]), Cvoid,
              (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ptr{Float64}),
              nr, ns, rmax, mesh)

        # Set solver type
        solver_i32 = Int32(solver)
        ccall((:py_set_solver_, libhprmat[]), Cvoid,
              (Ref{Int32},), solver_i32)

        new(nr, ns, rmax, solver, mesh)
    end
end

"""
    solve(solver, lval, qk, eta, cpot; solver_type=nothing)

Solve the R-matrix scattering problem.

# Arguments
- `solver::RMatrixSolver`: The solver instance
- `lval::Vector{Int32}`: Angular momentum for each channel
- `qk::Vector{Float64}`: Wave numbers (positive=open, negative=closed)
- `eta::Vector{Float64}`: Sommerfeld parameters
- `cpot::Array{ComplexF64,3}`: Coupling potentials (nr*ns, nch, nch)
- `solver_type::Int`: Override solver type (optional)

# Returns
- `cu::Matrix{ComplexF64}`: Collision (S) matrix
- `nopen::Int`: Number of open channels
"""
function solve(s::RMatrixSolver, lval::Vector{Int32}, qk::Vector{Float64},
               eta::Vector{Float64}, cpot::Array{ComplexF64,3};
               solver_type::Union{Nothing,Integer}=nothing)
    check_lib()

    nch = Int32(length(lval))
    nr = Int32(s.nr)
    ns = Int32(s.ns)
    rmax = s.rmax
    isolver = Int32(something(solver_type, s.solver))

    # Ensure Fortran memory layout (column-major)
    cpot_f = copy(cpot)

    # Output arrays
    cu = zeros(ComplexF64, nch, nch)
    nopen = Ref{Int32}(0)

    # Dummy arrays for wave function (not computed)
    cf_dummy = zeros(ComplexF64, 1, 1, 1)
    cpnl_dummy = zeros(ComplexF64, 1, 1, 1)
    nvc_dummy = Int32[1]
    twf = Int32(0)  # false
    ncp2 = Int32(0)

    # Call Fortran rmatrix
    ccall((:rmatrix_, libhprmat[]), Cvoid,
          (Ref{Int32},           # nch
           Ptr{Int32},           # lval
           Ptr{Float64},         # qk
           Ptr{Float64},         # eta
           Ref{Float64},         # rmax
           Ref{Int32},           # nr
           Ref{Int32},           # ns
           Ptr{ComplexF64},      # cpot
           Ptr{ComplexF64},      # cu
           Ref{Int32},           # ncp1
           Ref{Int32},           # ndim
           Ref{Int32},           # nopen
           Ref{Int32},           # twf (as integer for logical)
           Ptr{ComplexF64},      # cf
           Ref{Int32},           # nwf1
           Ref{Int32},           # nwf2
           Ref{Int32},           # nc
           Ptr{Int32},           # nvc
           Ref{Int32},           # ncp2
           Ptr{ComplexF64},      # cpnl
           Ref{Int32}),          # isolver
          nch, lval, qk, eta, rmax, nr, ns,
          cpot_f, cu, Int32(nr*ns), nch, nopen, twf,
          cf_dummy, Int32(1), Int32(1), Int32(1), nvc_dummy,
          ncp2, cpnl_dummy, isolver)

    return cu, Int(nopen[])
end

"""
    solve_with_wavefunction(solver, lval, qk, eta, cpot; entrance_channels=[1], solver_type=nothing)

Solve R-matrix problem and compute wave functions.

# Returns
- `cu::Matrix{ComplexF64}`: Collision matrix
- `cf::Array{ComplexF64,3}`: Wave functions (nr*ns, nch, n_entrance)
- `nopen::Int`: Number of open channels
"""
function solve_with_wavefunction(s::RMatrixSolver, lval::Vector{Int32},
                                  qk::Vector{Float64}, eta::Vector{Float64},
                                  cpot::Array{ComplexF64,3};
                                  entrance_channels::Vector{Int32}=Int32[1],
                                  solver_type::Union{Nothing,Integer}=nothing)
    check_lib()

    nch = Int32(length(lval))
    nr = Int32(s.nr)
    ns = Int32(s.ns)
    rmax = s.rmax
    isolver = Int32(something(solver_type, s.solver))
    nc_entrance = Int32(length(entrance_channels))

    cpot_f = copy(cpot)

    # Output arrays
    cu = zeros(ComplexF64, nch, nch)
    cf = zeros(ComplexF64, nr*ns, nch, nc_entrance)
    nopen = Ref{Int32}(0)

    cpnl_dummy = zeros(ComplexF64, 1, 1, 1)
    twf = Int32(1)  # true
    ncp2 = Int32(0)

    ccall((:rmatrix_, libhprmat[]), Cvoid,
          (Ref{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
           Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{ComplexF64},
           Ptr{ComplexF64}, Ref{Int32}, Ref{Int32}, Ref{Int32},
           Ref{Int32}, Ptr{ComplexF64}, Ref{Int32}, Ref{Int32},
           Ref{Int32}, Ptr{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ref{Int32}),
          nch, lval, qk, eta, rmax, nr, ns, cpot_f,
          cu, Int32(nr*ns), nch, nopen, twf, cf, Int32(nr*ns), nch,
          nc_entrance, entrance_channels, ncp2, cpnl_dummy, isolver)

    return cu, cf, Int(nopen[])
end

"""
    solve_nonlocal(solver, lval, qk, eta, cpot, cpnl; entrance_channels=[1], solver_type=nothing)

Solve R-matrix problem with non-local potential.
"""
function solve_nonlocal(s::RMatrixSolver, lval::Vector{Int32},
                         qk::Vector{Float64}, eta::Vector{Float64},
                         cpot::Array{ComplexF64,3}, cpnl::Array{ComplexF64,3};
                         entrance_channels::Vector{Int32}=Int32[1],
                         solver_type::Union{Nothing,Integer}=nothing)
    check_lib()

    nch = Int32(length(lval))
    nr = Int32(s.nr)
    ns = Int32(s.ns)
    rmax = s.rmax
    isolver = Int32(something(solver_type, s.solver))
    nc_entrance = Int32(length(entrance_channels))
    ncp2 = Int32(size(cpnl, 1))

    cpot_f = copy(cpot)
    cpnl_f = copy(cpnl)

    cu = zeros(ComplexF64, nch, nch)
    cf = zeros(ComplexF64, nr*ns, nch, nc_entrance)
    nopen = Ref{Int32}(0)
    twf = Int32(1)

    ccall((:rmatrix_, libhprmat[]), Cvoid,
          (Ref{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
           Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{ComplexF64},
           Ptr{ComplexF64}, Ref{Int32}, Ref{Int32}, Ref{Int32},
           Ref{Int32}, Ptr{ComplexF64}, Ref{Int32}, Ref{Int32},
           Ref{Int32}, Ptr{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ref{Int32}),
          nch, lval, qk, eta, rmax, nr, ns, cpot_f,
          cu, Int32(nr*ns), nch, nopen, twf, cf, Int32(nr*ns), nch,
          nc_entrance, entrance_channels, ncp2, cpnl_f, isolver)

    return cu, cf, Int(nopen[])
end

"""
    interpolate_wavefunction(solver, lval, qk, eta, cu, cf, nopen, channel, entrance, npoin; h=nothing)

Interpolate wave function onto uniform mesh.

# Returns
- `r::Vector{Float64}`: Radial coordinates
- `wf::Vector{ComplexF64}`: Wave function values
"""
function interpolate_wavefunction(s::RMatrixSolver, lval::Vector{Int32},
                                   qk::Vector{Float64}, eta::Vector{Float64},
                                   cu::Matrix{ComplexF64}, cf::Array{ComplexF64,3},
                                   nopen::Integer, channel::Integer,
                                   entrance::Integer, npoin::Integer;
                                   h::Union{Nothing,Real}=nothing)
    check_lib()

    nch = Int32(length(lval))
    nr = Int32(s.nr)
    ns = Int32(s.ns)
    nom = Int32(size(cf, 3))

    h_val = something(h, 2.0 * s.rmax / npoin)
    cwftab = zeros(ComplexF64, npoin)

    ccall((:wf_print_, libhprmat[]), Cvoid,
          (Ref{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
           Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{ComplexF64},
           Ref{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ref{Int32}, Ref{Int32},
           Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64},
           Ptr{ComplexF64}),
          nch, lval, qk, eta, s.rmax, nr, ns, cu, nch, Int32(nopen),
          cf, Int32(nr*ns), nch, s.mesh, Int32(channel), Int32(entrance),
          Int32(npoin), Float64(h_val), cwftab)

    r = collect(1:npoin) .* h_val
    return r, cwftab
end

"""
    set_solver!(solver, solver_type)

Set the default solver type.
"""
function set_solver!(s::RMatrixSolver, solver_type::Integer)
    check_lib()
    s.solver = solver_type
    solver_i32 = Int32(solver_type)
    ccall((:py_set_solver_, libhprmat[]), Cvoid, (Ref{Int32},), solver_i32)
    return s
end

"""
Get solver name.
"""
solver_name(s::RMatrixSolver) = get(SOLVER_NAMES, s.solver, "Unknown")

"""
List available solvers.
"""
available_solvers() = SOLVER_NAMES

end # module
