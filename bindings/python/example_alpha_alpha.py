#!/usr/bin/env python3
"""
Alpha-Alpha Scattering Example
==============================

This example demonstrates using HPRMAT for single-channel alpha-alpha
elastic scattering with the Ali-Bodmer potential.

Based on: Ex1 from the HPRMAT examples
Reference: Descouvemont CPC 200 (2016)
"""

import numpy as np
import sys

# Try to import the hprmat module
try:
    from hprmat import RMatrixSolver, SOLVER_DENSE, SOLVER_MIXED, SOLVER_WOODBURY
except ImportError:
    print("Error: hprmat module not found.")
    print("Build it first with: make")
    sys.exit(1)


def ali_bodmer_potential(r):
    """
    Ali-Bodmer alpha-alpha potential (deep version).

    V(r) = 500 * exp(-0.7^2 * r^2) - 130 * exp(-0.475^2 * r^2)  (MeV)
    """
    return 500.0 * np.exp(-0.7**2 * r**2) - 130.0 * np.exp(-0.475**2 * r**2)


def main():
    # Physical parameters
    mu = 2.0 * 4.0 / (2.0 + 4.0) * 938.0  # Reduced mass (MeV/c^2)
    # Actually for alpha-alpha: mu = 4 * 938 / 2 = 1876 MeV
    mu = 4.0 * 938.0 / 2.0
    hbarc = 197.3  # MeV.fm

    # R-matrix parameters
    nr = 30  # Lagrange functions
    ns = 1   # Single interval
    rmax = 10.0  # Channel radius (fm)

    # Create solver
    print("Initializing R-matrix solver...")
    solver = RMatrixSolver(nr, ns, rmax, solver=SOLVER_DENSE)
    print(f"  Mesh points: {len(solver.mesh)}")
    print(f"  Solver: {solver.solver_name}")

    # Channel parameters (single channel, L=0)
    nch = 1
    lval = np.array([0], dtype=np.int32)  # s-wave
    eta = np.array([0.0], dtype=np.float64)  # No Coulomb (uncharged)

    # Build potential matrix on Lagrange mesh
    cpot = np.zeros((nr * ns, nch, nch), dtype=np.complex128, order='F')
    for ir, r in enumerate(solver.mesh):
        cpot[ir, 0, 0] = 2.0 * mu / hbarc**2 * ali_bodmer_potential(r)

    # Energy scan
    energies = [1.0, 2.0, 3.0, 4.0, 5.0]  # MeV (CM)

    print("\nAlpha-Alpha Elastic Scattering (Ali-Bodmer Potential)")
    print("=" * 60)
    print(f"{'E_cm (MeV)':>12} {'k (fm^-1)':>12} {'Re(S)':>12} {'Im(S)':>12} {'delta (deg)':>12}")
    print("-" * 60)

    for E in energies:
        # Wave number
        k = np.sqrt(2.0 * mu * E) / hbarc  # fm^-1
        qk = np.array([k], dtype=np.float64)

        # Solve
        S, nopen = solver.solve(lval, qk, eta, cpot)

        # Extract phase shift from S-matrix
        # S = exp(2i*delta) for elastic scattering
        S_11 = S[0, 0]
        delta = np.angle(S_11) / 2.0 * 180.0 / np.pi

        print(f"{E:12.2f} {k:12.4f} {S_11.real:12.6f} {S_11.imag:12.6f} {delta:12.2f}")

    # Compare solvers
    print("\n\nSolver Comparison at E = 4 MeV")
    print("=" * 60)

    E = 4.0
    k = np.sqrt(2.0 * mu * E) / hbarc
    qk = np.array([k], dtype=np.float64)

    solvers = [
        (SOLVER_DENSE, "Dense LAPACK"),
        (SOLVER_MIXED, "Mixed Precision"),
        (SOLVER_WOODBURY, "Woodbury"),
    ]

    results = []
    for stype, name in solvers:
        S, nopen = solver.solve(lval, qk, eta, cpot, solver=stype)
        S_11 = S[0, 0]
        results.append((name, S_11))
        print(f"{name:20s}: S = {S_11.real:12.8f} + {S_11.imag:12.8f}i")

    # Show differences
    ref = results[0][1]
    print("\nDifference from Dense LAPACK:")
    for name, S_11 in results[1:]:
        diff = abs(S_11 - ref)
        print(f"  {name:20s}: |dS| = {diff:.2e}")

    print("\nDone!")


if __name__ == "__main__":
    main()
