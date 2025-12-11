#!/usr/bin/env julia
"""
Test Ex1: alpha-208Pb Optical Potential (Julia)
"""

include("HPRMAT.jl")
using .HPRMAT
using Printf

function main()
    println("=" ^ 60)
    println("HPRMAT Julia Test: Ex1 alpha-208Pb Optical Potential")
    println("=" ^ 60)

    # Physical constants
    a1, a2 = 208.0, 4.0
    z1, z2 = 82.0, 2.0
    rmu = a1 * a2 / (a1 + a2)
    ze = z1 * z2 * 1.44
    hm = 20.736 / rmu
    eta0 = ze / (2 * hm)

    # R-matrix parameters
    l, nr, ns, rmax = 20, 60, 1, 14.0
    ne, e0, estep = 5, 10.0, 10.0

    # Woods-Saxon parameters
    an = 0.5803
    rn = 1.1132 * (a1^(1/3) + a2^(1/3))

    println("\nParameters: L=$l, nr=$nr, ns=$ns, rmax=$rmax fm")

    # Initialize
    zrma = rmat_init(nr, ns, rmax)
    ntot = nr * ns

    # Build potential
    nch = 1
    cpot = zeros(ComplexF64, ntot, nch, nch)
    for i in 1:ntot
        r = zrma[i]
        cvn = complex(-100.0, -10.0) / (1.0 + exp((r - rn) / an))
        cpot[i, 1, 1] = (cvn + ze / r) / hm
    end

    lval = Int32[l]

    println("-" ^ 60)
    @printf("%10s %14s %14s %10s\n", "E (MeV)", "Re(S)", "Im(S)", "|S|")
    println("-" ^ 60)

    # Test solvers
    for (solver, name) in [(1,"Dense"), (2,"Mixed"), (3,"Woodbury")]
        println("\n--- Solver $solver: $name ---")
        for ie in 0:(ne-1)
            ecm = e0 + ie * estep
            qk = [sqrt(ecm / hm)]
            eta = [eta0 / qk[1]]

            cu, nopen = rmatrix(nch, lval, qk, eta, rmax, nr, ns, cpot, solver)
            S = cu[1,1]
            @printf("%10.3f %14.6e %14.6e %10.6f\n", ecm, real(S), imag(S), abs(S))
        end
    end

    println("-" ^ 60)
    println("\nTest completed!")
end

main()
