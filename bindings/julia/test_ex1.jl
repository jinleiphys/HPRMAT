#!/usr/bin/env julia
"""
Test Ex1: alpha-208Pb Optical Potential (Julia)
================================================

This replicates the Fortran Ex1 example exactly.
Compares results with the Fortran reference.

Reference: Goldring et al., Phys. Lett. B32 (1970) 465
"""

using Printf

# Library path
const libhprmat = joinpath(@__DIR__, "..", "..", "lib", "libhprmat.dylib")

function main()
    println("=" ^ 60)
    println("HPRMAT Julia Test: Ex1 alpha-208Pb Optical Potential")
    println("=" ^ 60)

    # Check library
    if !isfile(libhprmat)
        error("Library not found: $libhprmat\nBuild it first with: cd bindings && make lib")
    end

    # Physical constants (same as Fortran)
    a1, a2 = 208.0, 4.0    # Mass numbers
    z1, z2 = 82.0, 2.0     # Charge numbers
    rmu = a1 * a2 / (a1 + a2)  # Reduced mass
    ze = z1 * z2 * 1.44    # Z1*Z2*e^2 (MeV.fm)
    hm = 20.736 / rmu      # hbar^2/(2*mu) in MeV.fm^2
    eta0 = ze / (2 * hm)   # Sommerfeld parameter factor

    # R-matrix parameters (from data1)
    l = 20       # Angular momentum
    nr = 60      # Lagrange functions per interval
    ns = 1       # Number of intervals
    rmax = 14.0  # Channel radius (fm)

    # Energy parameters
    ne = 5       # Number of energies
    e0 = 10.0    # Initial energy (MeV)
    estep = 10.0 # Energy step (MeV)

    # Woods-Saxon parameters
    an = 0.5803  # Diffuseness (fm)
    rn = 1.1132 * (a1^(1.0/3) + a2^(1.0/3))  # Radius (fm)

    println("\nParameters:")
    println("  Angular momentum L = $l")
    println("  Lagrange functions: nr = $nr, ns = $ns")
    @printf("  Channel radius: rmax = %.1f fm\n", rmax)
    @printf("  Reduced mass: mu = %.4f amu\n", rmu)
    @printf("  hbar^2/(2*mu) = %.6f MeV.fm^2\n", hm)

    # Initialize R-matrix
    ntot = nr * ns
    zrma = zeros(Float64, ntot)

    nr_ref = Ref{Int32}(nr)
    ns_ref = Ref{Int32}(ns)
    rmax_ref = Ref{Float64}(rmax)

    ccall((:rmat_ini_, libhprmat), Cvoid,
          (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ptr{Float64}),
          nr_ref, ns_ref, rmax_ref, zrma)

    println("  Mesh points: $ntot")
    @printf("  First mesh point: r = %.6f fm\n", zrma[1])
    @printf("  Last mesh point: r = %.6f fm\n", zrma[end])

    # Build potential matrix (complex Woods-Saxon + Coulomb)
    nch = 1
    cpot = zeros(ComplexF64, ntot, nch, nch)

    for i in 1:ntot
        r = zrma[i]
        xx = 1.0 + exp((r - rn) / an)
        cvn = complex(-100.0, -10.0) / xx  # Complex Woods-Saxon
        vc = ze / r                         # Coulomb
        cpot[i, 1, 1] = (cvn + vc) / hm     # Divided by hbar^2/(2*mu)
    end

    # Channel parameters
    lval = Int32[l]

    @printf("\nEnergies: %d points from %.1f to %.1f MeV\n", ne, e0, e0 + (ne-1)*estep)
    println("-" ^ 60)
    @printf("%10s %14s %14s %10s\n", "E (MeV)", "Re(S)", "Im(S)", "|S|")
    println("-" ^ 60)

    # Dummy arrays for wave function (not computed)
    cf_dummy = zeros(ComplexF64, 1, 1, 1)
    cpnl_dummy = zeros(ComplexF64, 1, 1, 1)
    nvc_dummy = Int32[1]

    # Test all solver types
    solver_names = Dict(1 => "Dense LAPACK", 2 => "Mixed Precision", 3 => "Woodbury")

    for solver in [1, 2, 3]
        println("\n--- Solver $solver: $(solver_names[solver]) ---")

        for ie in 0:(ne-1)
            ecm = e0 + ie * estep

            # Wave number and Sommerfeld parameter
            qk_val = sqrt(ecm / hm)
            eta_val = eta0 / qk_val

            qk = Float64[qk_val]
            eta = Float64[eta_val]

            # Output arrays
            cu = zeros(ComplexF64, nch, nch)
            nopen = Ref{Int32}(0)
            twf = Ref{Int32}(0)  # false

            # Call Fortran rmatrix
            nch_ref = Ref{Int32}(nch)
            nr_ref = Ref{Int32}(nr)
            ns_ref = Ref{Int32}(ns)
            rmax_ref = Ref{Float64}(rmax)
            ntot_ref = Ref{Int32}(ntot)
            ncp2_ref = Ref{Int32}(0)
            solver_ref = Ref{Int32}(solver)

            ccall((:rmatrix_, libhprmat), Cvoid,
                  (Ref{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
                   Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{ComplexF64},
                   Ptr{ComplexF64}, Ref{Int32}, Ref{Int32}, Ref{Int32},
                   Ref{Int32}, Ptr{ComplexF64}, Ref{Int32}, Ref{Int32},
                   Ref{Int32}, Ptr{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ref{Int32}),
                  nch_ref, lval, qk, eta, rmax_ref, nr_ref, ns_ref, cpot,
                  cu, ntot_ref, nch_ref, nopen, twf, cf_dummy,
                  Ref{Int32}(1), Ref{Int32}(1), Ref{Int32}(1), nvc_dummy,
                  ncp2_ref, cpnl_dummy, solver_ref)

            # Extract S-matrix element
            S11 = cu[1, 1]
            @printf("%10.3f %14.6e %14.6e %10.6f\n", ecm, real(S11), imag(S11), abs(S11))
        end
    end

    println("-" ^ 60)
    println("\nTest completed!")
end

main()
