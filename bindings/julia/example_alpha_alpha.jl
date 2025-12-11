#!/usr/bin/env julia
"""
Alpha-Alpha Scattering Example (Julia)
======================================

This example demonstrates using HPRMAT.jl for single-channel alpha-alpha
elastic scattering with the Ali-Bodmer potential.

Based on: Ex1 from the HPRMAT examples
Reference: Descouvemont CPC 200 (2016)

Run with: julia example_alpha_alpha.jl
"""

# Include the HPRMAT module
include("HPRMAT.jl")
using .HPRMAT
using Printf

"""
Ali-Bodmer alpha-alpha potential (deep version).

V(r) = 500 * exp(-0.7^2 * r^2) - 130 * exp(-0.475^2 * r^2)  (MeV)
"""
function ali_bodmer_potential(r)
    return 500.0 * exp(-0.7^2 * r^2) - 130.0 * exp(-0.475^2 * r^2)
end

function main()
    # Physical parameters
    mu = 4.0 * 938.0 / 2.0  # Reduced mass for alpha-alpha (MeV/c^2)
    hbarc = 197.3  # MeV.fm

    # R-matrix parameters
    nr = 30  # Lagrange functions
    ns = 1   # Single interval
    rmax = 10.0  # Channel radius (fm)

    # Create solver
    println("Initializing R-matrix solver...")
    solver = RMatrixSolver(nr, ns, rmax; solver=SOLVER_DENSE)
    println("  Mesh points: $(length(solver.mesh))")
    println("  Solver: $(solver_name(solver))")

    # Channel parameters (single channel, L=0)
    nch = 1
    lval = Int32[0]  # s-wave
    eta = [0.0]      # No Coulomb

    # Build potential matrix on Lagrange mesh
    cpot = zeros(ComplexF64, nr * ns, nch, nch)
    for (ir, r) in enumerate(solver.mesh)
        cpot[ir, 1, 1] = 2.0 * mu / hbarc^2 * ali_bodmer_potential(r)
    end

    # Energy scan
    energies = [1.0, 2.0, 3.0, 4.0, 5.0]  # MeV (CM)

    println("\nAlpha-Alpha Elastic Scattering (Ali-Bodmer Potential)")
    println("=" ^ 60)
    @printf("%12s %12s %12s %12s %12s\n", "E_cm (MeV)", "k (fm^-1)", "Re(S)", "Im(S)", "delta (deg)")
    println("-" ^ 60)

    for E in energies
        # Wave number
        k = sqrt(2.0 * mu * E) / hbarc  # fm^-1
        qk = [k]

        # Solve
        S, nopen = solve(solver, lval, qk, eta, cpot)

        # Extract phase shift from S-matrix
        S_11 = S[1, 1]
        delta = angle(S_11) / 2.0 * 180.0 / π

        @printf("%12.2f %12.4f %12.6f %12.6f %12.2f\n", E, k, real(S_11), imag(S_11), delta)
    end

    # Compare solvers
    println("\n\nSolver Comparison at E = 4 MeV")
    println("=" ^ 60)

    E = 4.0
    k = sqrt(2.0 * mu * E) / hbarc
    qk = [k]

    solvers_list = [
        (SOLVER_DENSE, "Dense LAPACK"),
        (SOLVER_MIXED, "Mixed Precision"),
        (SOLVER_WOODBURY, "Woodbury"),
    ]

    results = Tuple{String, ComplexF64}[]
    for (stype, name) in solvers_list
        S, nopen = solve(solver, lval, qk, eta, cpot; solver_type=stype)
        S_11 = S[1, 1]
        push!(results, (name, S_11))
        @printf("%20s: S = %12.8f + %12.8fi\n", name, real(S_11), imag(S_11))
    end

    # Show differences
    ref = results[1][2]
    println("\nDifference from Dense LAPACK:")
    for (name, S_11) in results[2:end]
        diff = abs(S_11 - ref)
        @printf("  %20s: |dS| = %.2e\n", name, diff)
    end

    println("\nDone!")
end

# Run main function
main()
