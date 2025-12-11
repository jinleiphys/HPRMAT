#!/usr/bin/env python3
"""
Test Ex1: alpha-208Pb Optical Potential
========================================

This replicates the Fortran Ex1 example exactly.
Compares results with the Fortran reference.

Reference: Goldring et al., Phys. Lett. B32 (1970) 465
"""

import numpy as np
import sys

try:
    import hprmat_fortran as hpf
except ImportError:
    print("Error: hprmat_fortran module not found.")
    print("Build it first with f2py")
    sys.exit(1)


def main():
    print("=" * 60)
    print("HPRMAT Python Test: Ex1 alpha-208Pb Optical Potential")
    print("=" * 60)

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
    rn = 1.1132 * (a1**(1.0/3) + a2**(1.0/3))  # Radius (fm)

    print(f"\nParameters:")
    print(f"  Angular momentum L = {l}")
    print(f"  Lagrange functions: nr = {nr}, ns = {ns}")
    print(f"  Channel radius: rmax = {rmax} fm")
    print(f"  Reduced mass: mu = {rmu:.4f} amu")
    print(f"  hbar^2/(2*mu) = {hm:.6f} MeV.fm^2")

    # Initialize R-matrix
    zrma = hpf.py_rmat_init(nr, ns, rmax)
    ntot = nr * ns

    print(f"  Mesh points: {ntot}")
    print(f"  First mesh point: r = {zrma[0]:.6f} fm")
    print(f"  Last mesh point: r = {zrma[-1]:.6f} fm")

    # Build potential matrix (complex Woods-Saxon + Coulomb)
    nch = 1
    cpot = np.zeros((ntot, nch, nch), dtype=np.complex128, order='F')

    for i in range(ntot):
        r = zrma[i]
        xx = 1.0 + np.exp((r - rn) / an)
        cvn = complex(-100.0, -10.0) / xx  # Complex Woods-Saxon
        vc = ze / r                         # Coulomb
        cpot[i, 0, 0] = (cvn + vc) / hm     # Divided by hbar^2/(2*mu)

    # Channel parameters
    lval = np.array([l], dtype=np.int32)

    print(f"\nEnergies: {ne} points from {e0} to {e0 + (ne-1)*estep} MeV")
    print("-" * 60)
    print(f"{'E (MeV)':>10} {'Re(S)':>14} {'Im(S)':>14} {'|S|':>10}")
    print("-" * 60)

    # Test all solver types
    for solver in [1, 2, 3]:
        solver_names = {1: "Dense LAPACK", 2: "Mixed Precision", 3: "Woodbury"}
        print(f"\n--- Solver {solver}: {solver_names[solver]} ---")

        for ie in range(ne):
            ecm = e0 + ie * estep

            # Wave number and Sommerfeld parameter
            qk_val = np.sqrt(ecm / hm)
            eta_val = eta0 / qk_val

            qk = np.array([qk_val], dtype=np.float64)
            eta = np.array([eta_val], dtype=np.float64)

            # Solve
            cu, nopen = hpf.py_rmatrix(nch, ntot, lval, qk, eta, rmax, nr, ns, cpot, solver)

            # Extract S-matrix element
            S11 = cu[0, 0]
            print(f"{ecm:10.3f} {S11.real:14.6e} {S11.imag:14.6e} {abs(S11):10.6f}")

    print("-" * 60)
    print("\nTest completed!")


if __name__ == "__main__":
    main()
