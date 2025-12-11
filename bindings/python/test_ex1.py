#!/usr/bin/env python3
"""
Test Ex1: alpha-208Pb Optical Potential
========================================

This example demonstrates how to use HPRMAT Python bindings to calculate
elastic scattering of alpha particles on 208Pb using a complex Woods-Saxon
optical potential.

Physical system:
    - Projectile: alpha (4He), A=4, Z=2
    - Target: 208Pb, A=208, Z=82
    - Potential: Complex Woods-Saxon + Coulomb

Reference: Goldring et al., Phys. Lett. B32 (1970) 465
"""

import numpy as np
import sys

# =============================================================================
# Import the f2py-generated Fortran module
# Build it first with: make (in the python directory)
# =============================================================================
try:
    import hprmat_fortran as hpf
except ImportError:
    print("Error: hprmat_fortran module not found.")
    print("Build it first with: make")
    sys.exit(1)


def main():
    print("=" * 60)
    print("HPRMAT Python Test: Ex1 alpha-208Pb Optical Potential")
    print("=" * 60)

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
    rn = 1.1132 * (a1**(1.0/3) + a2**(1.0/3))  # Radius parameter in fm

    print(f"\nParameters:")
    print(f"  Angular momentum L = {l}")
    print(f"  Lagrange mesh: nr={nr}, ns={ns}, total points={nr*ns}")
    print(f"  Channel radius: rmax = {rmax} fm")

    # =========================================================================
    # Step 1: Initialize R-matrix calculation
    # This computes Lagrange mesh points and kinetic energy matrix elements
    # Returns: zrma = array of radial mesh points (size: nr*ns)
    # =========================================================================
    zrma = hpf.py_rmat_init(nr, ns, rmax)
    ntot = nr * ns

    print(f"  Mesh range: {zrma[0]:.4f} to {zrma[-1]:.4f} fm")

    # =========================================================================
    # Step 2: Build the coupling potential matrix
    #
    # cpot[ir, i, j] = V_ij(r_ir) / (hbar^2/2mu)
    #
    # For single channel (nch=1): only cpot[:, 0, 0] is used
    # The potential must be divided by hbar^2/(2*mu) for HPRMAT convention
    #
    # Here we use complex Woods-Saxon + Coulomb:
    #   V(r) = V_WS(r) + V_C(r)
    #   V_WS = -(V0 + i*W0) / (1 + exp((r-R)/a))  with V0=100, W0=10 MeV
    #   V_C  = Z1*Z2*e^2 / r
    # =========================================================================
    nch = 1  # Number of channels (single channel for elastic scattering)
    cpot = np.zeros((ntot, nch, nch), dtype=np.complex128, order='F')

    for i in range(ntot):
        r = zrma[i]
        # Complex Woods-Saxon potential (absorptive)
        xx = 1.0 + np.exp((r - rn) / an)
        cvn = complex(-100.0, -10.0) / xx  # V0=100 MeV, W0=10 MeV
        # Point Coulomb potential
        vc = ze / r
        # Store potential divided by hbar^2/(2*mu)
        cpot[i, 0, 0] = (cvn + vc) / hm

    # =========================================================================
    # Step 3: Set up channel quantum numbers
    #
    # lval[i] = orbital angular momentum for channel i
    # For single channel elastic scattering, lval = [L]
    # =========================================================================
    lval = np.array([l], dtype=np.int32)

    print(f"\nEnergy scan: {ne} points from {e0} to {e0 + (ne-1)*estep} MeV")
    print("-" * 60)
    print(f"{'E (MeV)':>10} {'Re(S)':>14} {'Im(S)':>14} {'|S|':>10}")
    print("-" * 60)

    # =========================================================================
    # Step 4: Loop over energies and solver types
    #
    # For each energy E:
    #   - Calculate wave number: k = sqrt(2*mu*E) / hbar = sqrt(E/hm)
    #   - Calculate Sommerfeld parameter: eta = Z1*Z2*e^2*mu / (hbar^2*k)
    #   - Call py_rmatrix to get the collision (S) matrix
    #
    # Solver types:
    #   1 = Dense LAPACK (ZGESV) - reference, most accurate
    #   2 = Mixed precision - faster, slightly less accurate
    #   3 = Woodbury - CPU optimized, good accuracy
    # =========================================================================
    for solver in [1, 2, 3]:
        solver_names = {1: "Dense LAPACK", 2: "Mixed Precision", 3: "Woodbury"}
        print(f"\n--- Solver {solver}: {solver_names[solver]} ---")

        for ie in range(ne):
            ecm = e0 + ie * estep  # Center-of-mass energy in MeV

            # Wave number k in fm^-1
            qk_val = np.sqrt(ecm / hm)
            # Sommerfeld parameter (dimensionless)
            eta_val = eta0 / qk_val

            # Arrays for HPRMAT (size = number of channels)
            qk = np.array([qk_val], dtype=np.float64)
            eta = np.array([eta_val], dtype=np.float64)

            # =================================================================
            # Call the R-matrix solver
            #
            # Arguments:
            #   nch    - number of channels
            #   ntot   - total mesh points (nr*ns)
            #   lval   - angular momentum array
            #   qk     - wave number array (positive=open, negative=closed)
            #   eta    - Sommerfeld parameter array
            #   rmax   - channel radius
            #   nr, ns - mesh parameters
            #   cpot   - coupling potential matrix
            #   solver - solver type (1, 2, or 3)
            #
            # Returns:
            #   cu    - collision (S) matrix, size (nch, nch)
            #   nopen - number of open channels
            # =================================================================
            cu, nopen = hpf.py_rmatrix(nch, ntot, lval, qk, eta,
                                        rmax, nr, ns, cpot, solver)

            # =================================================================
            # Extract results
            #
            # For elastic scattering, S-matrix element S_11 gives:
            #   |S| = 1 for no absorption
            #   |S| < 1 indicates absorption (flux loss)
            #   Phase shift: delta = arg(S) / 2
            # =================================================================
            S11 = cu[0, 0]
            print(f"{ecm:10.3f} {S11.real:14.6e} {S11.imag:14.6e} {abs(S11):10.6f}")

    print("-" * 60)
    print("\nTest completed!")


if __name__ == "__main__":
    main()
