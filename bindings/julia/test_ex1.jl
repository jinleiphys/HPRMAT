#!/usr/bin/env julia
"""
Test Ex1: alpha-208Pb Optical Potential (Julia)
================================================

This example demonstrates how to use HPRMAT Julia bindings to calculate
elastic scattering of alpha particles on 208Pb using a complex Woods-Saxon
optical potential.

Physical system:
    - Projectile: alpha (4He), A=4, Z=2
    - Target: 208Pb, A=208, Z=82
    - Potential: Complex Woods-Saxon + Coulomb

Reference: Goldring et al., Phys. Lett. B32 (1970) 465
"""

# =============================================================================
# Load the HPRMAT module
# The module provides: rmat_init(), rmatrix(), and solver constants
# =============================================================================
include("src/HPRMAT.jl")
using .HPRMAT
using Printf

function main()
    println("=" ^ 60)
    println("HPRMAT Julia Test: Ex1 alpha-208Pb Optical Potential")
    println("=" ^ 60)

    # =========================================================================
    # Physical constants
    # =========================================================================
    a1, a2 = 208.0, 4.0    # Mass numbers (target, projectile)
    z1, z2 = 82.0, 2.0     # Charge numbers (target, projectile)
    rmu = a1 * a2 / (a1 + a2)  # Reduced mass in amu
    ze = z1 * z2 * 1.44    # Z1*Z2*e^2 in MeV.fm (Coulomb constant)
    hm = 20.736 / rmu      # hbar^2/(2*mu) in MeV.fm^2
    eta0 = ze / (2 * hm)   # Sommerfeld parameter prefactor: eta = eta0/k

    # =========================================================================
    # R-matrix parameters
    # =========================================================================
    l = 20       # Total angular momentum (partial wave)
    nr = 60      # Number of Lagrange-Legendre mesh points per interval
    ns = 1       # Number of radial intervals (usually 1)
    rmax = 14.0  # Channel radius in fm (boundary of internal region)

    # Energy scan parameters
    ne = 5       # Number of energy points
    e0 = 10.0    # Starting energy in MeV (center-of-mass)
    estep = 10.0 # Energy step in MeV

    # Woods-Saxon potential parameters
    an = 0.5803  # Diffuseness in fm
    rn = 1.1132 * (a1^(1/3) + a2^(1/3))  # Radius parameter in fm

    println("\nParameters: L=$l, nr=$nr, ns=$ns, rmax=$rmax fm")

    # =========================================================================
    # Step 1: Initialize R-matrix calculation
    #
    # rmat_init(nr, ns, rmax) computes:
    #   - Lagrange-Legendre mesh points
    #   - Kinetic energy matrix elements
    #   - Bloch operator matrix elements
    #
    # Returns: zrma = array of radial mesh points (size: nr*ns)
    # =========================================================================
    zrma = rmat_init(nr, ns, rmax)
    ntot = nr * ns

    # =========================================================================
    # Step 2: Build the coupling potential matrix
    #
    # cpot[ir, i, j] = V_ij(r_ir) / (hbar^2/2mu)
    #
    # For single channel (nch=1): only cpot[:, 1, 1] is used (Julia is 1-indexed)
    # The potential must be divided by hbar^2/(2*mu) for HPRMAT convention
    #
    # Here we use complex Woods-Saxon + Coulomb:
    #   V(r) = V_WS(r) + V_C(r)
    #   V_WS = -(V0 + i*W0) / (1 + exp((r-R)/a))  with V0=100, W0=10 MeV
    #   V_C  = Z1*Z2*e^2 / r
    # =========================================================================
    nch = 1  # Number of channels (single channel for elastic scattering)
    cpot = zeros(ComplexF64, ntot, nch, nch)

    for i in 1:ntot
        r = zrma[i]
        # Complex Woods-Saxon potential (absorptive)
        cvn = complex(-100.0, -10.0) / (1.0 + exp((r - rn) / an))
        # Point Coulomb + Woods-Saxon, divided by hbar^2/(2*mu)
        cpot[i, 1, 1] = (cvn + ze / r) / hm
    end

    # =========================================================================
    # Step 3: Set up channel quantum numbers
    #
    # lval[i] = orbital angular momentum for channel i
    # Must be Int32 for Fortran compatibility
    # =========================================================================
    lval = Int32[l]

    println("-" ^ 60)
    @printf("%10s %14s %14s %10s\n", "E (MeV)", "Re(S)", "Im(S)", "|S|")
    println("-" ^ 60)

    # =========================================================================
    # Step 4: Loop over energies and solver types
    #
    # For each energy E:
    #   - Calculate wave number: k = sqrt(2*mu*E) / hbar = sqrt(E/hm)
    #   - Calculate Sommerfeld parameter: eta = Z1*Z2*e^2*mu / (hbar^2*k)
    #   - Call rmatrix() to get the collision (S) matrix
    #
    # Solver types (from HPRMAT module):
    #   SOLVER_DENSE (1)    - Dense LAPACK, reference accuracy
    #   SOLVER_MIXED (2)    - Mixed precision, faster
    #   SOLVER_WOODBURY (3) - CPU optimized
    #   SOLVER_GPU (4)      - GPU accelerated (if available)
    # =========================================================================
    for (solver, name) in [(SOLVER_DENSE, "Dense LAPACK"),
                           (SOLVER_MIXED, "Mixed Precision"),
                           (SOLVER_WOODBURY, "Woodbury")]
        println("\n--- Solver $solver: $name ---")

        for ie in 0:(ne-1)
            ecm = e0 + ie * estep  # Center-of-mass energy in MeV

            # Wave number k in fm^-1
            qk_val = sqrt(ecm / hm)
            # Sommerfeld parameter (dimensionless)
            eta_val = eta0 / qk_val

            # Arrays for HPRMAT (size = number of channels)
            qk = [qk_val]
            eta = [eta_val]

            # =================================================================
            # Call the R-matrix solver
            #
            # rmatrix(nch, lval, qk, eta, rmax, nr, ns, cpot, solver)
            #
            # Arguments:
            #   nch    - number of channels
            #   lval   - angular momentum array (Int32)
            #   qk     - wave number array (positive=open, negative=closed)
            #   eta    - Sommerfeld parameter array
            #   rmax   - channel radius
            #   nr, ns - mesh parameters
            #   cpot   - coupling potential matrix (ComplexF64)
            #   solver - solver type (1, 2, 3, or 4)
            #
            # Returns:
            #   cu    - collision (S) matrix, size (nch, nch)
            #   nopen - number of open channels
            # =================================================================
            cu, nopen = rmatrix(nch, lval, qk, eta, rmax, nr, ns, cpot, solver)

            # =================================================================
            # Extract results
            #
            # For elastic scattering, S-matrix element S_11 gives:
            #   |S| = 1 for no absorption
            #   |S| < 1 indicates absorption (flux loss)
            #   Phase shift: delta = angle(S) / 2
            # =================================================================
            S = cu[1, 1]  # Julia is 1-indexed
            @printf("%10.3f %14.6e %14.6e %10.6f\n", ecm, real(S), imag(S), abs(S))
        end
    end

    println("-" ^ 60)
    println("\nTest completed!")
end

main()
