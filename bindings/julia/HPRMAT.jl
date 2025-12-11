"""
HPRMAT.jl - High-Performance R-Matrix Solver for Julia

Simple interface to the HPRMAT Fortran library.

# Example
```julia
include("HPRMAT.jl")
using .HPRMAT

# Initialize
zrma = rmat_init(60, 1, 14.0)

# Build potential and solve
cu, nopen = rmatrix(nch, lval, qk, eta, rmax, nr, ns, cpot, solver)
```
"""
module HPRMAT

export rmat_init, rmatrix, set_solver
export SOLVER_DENSE, SOLVER_MIXED, SOLVER_WOODBURY, SOLVER_GPU

# Solver constants
const SOLVER_DENSE = 1
const SOLVER_MIXED = 2
const SOLVER_WOODBURY = 3
const SOLVER_GPU = 4

# Find library
function find_lib()
    paths = [
        joinpath(@__DIR__, "..", "..", "lib", "libhprmat.dylib"),
        joinpath(@__DIR__, "..", "..", "lib", "libhprmat.so"),
        joinpath(@__DIR__, "..", "lib", "libhprmat.dylib"),
        joinpath(@__DIR__, "..", "lib", "libhprmat.so"),
    ]
    for p in paths
        isfile(p) && return p
    end
    error("libhprmat not found. Build with: cd bindings && make lib")
end

const libhprmat = find_lib()

"""
    rmat_init(nr, ns, rmax) -> zrma

Initialize R-matrix calculation. Returns mesh points.
"""
function rmat_init(nr::Int, ns::Int, rmax::Float64)
    zrma = zeros(Float64, nr * ns)
    ccall((:rmat_ini_, libhprmat), Cvoid,
          (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ptr{Float64}),
          Int32(nr), Int32(ns), rmax, zrma)
    return zrma
end

"""
    rmatrix(nch, lval, qk, eta, rmax, nr, ns, cpot, solver) -> (cu, nopen)

Solve R-matrix problem. Returns collision matrix and number of open channels.
"""
function rmatrix(nch::Int, lval::Vector{Int32}, qk::Vector{Float64},
                 eta::Vector{Float64}, rmax::Float64, nr::Int, ns::Int,
                 cpot::Array{ComplexF64,3}, solver::Int=SOLVER_DENSE)

    ntot = nr * ns
    cu = zeros(ComplexF64, nch, nch)
    nopen = Ref{Int32}(0)

    # Dummy arrays
    cf = zeros(ComplexF64, 1, 1, 1)
    cpnl = zeros(ComplexF64, 1, 1, 1)
    nvc = Int32[1]

    ccall((:rmatrix_, libhprmat), Cvoid,
          (Ref{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
           Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{ComplexF64},
           Ptr{ComplexF64}, Ref{Int32}, Ref{Int32}, Ref{Int32},
           Ref{Int32}, Ptr{ComplexF64}, Ref{Int32}, Ref{Int32},
           Ref{Int32}, Ptr{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ref{Int32}),
          Int32(nch), lval, qk, eta, rmax, Int32(nr), Int32(ns), cpot,
          cu, Int32(ntot), Int32(nch), nopen, Int32(0), cf,
          Int32(1), Int32(1), Int32(1), nvc, Int32(0), cpnl, Int32(solver))

    return cu, Int(nopen[])
end

"""
    set_solver(solver)

Set default solver type (1-4).
"""
function set_solver(solver::Int)
    ccall((:py_set_solver_, libhprmat), Cvoid, (Ref{Int32},), Int32(solver))
end

end # module
