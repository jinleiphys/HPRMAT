/*
 * Test Ex1: alpha-208Pb Optical Potential (C)
 * ============================================
 *
 * This example demonstrates how to use the HPRMAT C interface to calculate
 * elastic scattering of alpha particles on 208Pb using a complex Woods-Saxon
 * optical potential.
 *
 * Physical system:
 *     - Projectile: alpha (4He), A=4, Z=2
 *     - Target: 208Pb, A=208, Z=82
 *     - Potential: Complex Woods-Saxon + Coulomb
 *
 * Reference: Goldring et al., Phys. Lett. B32 (1970) 465
 *
 * Build:
 *   gcc -o test_ex1 test_ex1.c -L../../lib -lhprmat -lm
 *
 * Run:
 *   DYLD_LIBRARY_PATH=../../lib ./test_ex1   (macOS)
 *   LD_LIBRARY_PATH=../../lib ./test_ex1     (Linux)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "hprmat.h"

int main() {
    printf("============================================================\n");
    printf("HPRMAT C Test: Ex1 alpha-208Pb Optical Potential\n");
    printf("============================================================\n");

    /* ========================================================================
     * Physical constants
     * ======================================================================== */
    double a1 = 208.0, a2 = 4.0;    /* Mass numbers (target, projectile) */
    double z1 = 82.0, z2 = 2.0;     /* Charge numbers (target, projectile) */
    double rmu = a1 * a2 / (a1 + a2);  /* Reduced mass in amu */
    double ze = z1 * z2 * 1.44;     /* Z1*Z2*e^2 in MeV.fm (Coulomb constant) */
    double hm = 20.736 / rmu;       /* hbar^2/(2*mu) in MeV.fm^2 */
    double eta0 = ze / (2 * hm);    /* Sommerfeld parameter prefactor */

    /* ========================================================================
     * R-matrix parameters
     * ======================================================================== */
    int l = 20;         /* Total angular momentum (partial wave) */
    int nr = 60;        /* Number of Lagrange-Legendre mesh points per interval */
    int ns = 1;         /* Number of radial intervals (usually 1) */
    double rmax = 14.0; /* Channel radius in fm (boundary of internal region) */

    /* Energy scan parameters */
    int ne = 5;           /* Number of energy points */
    double e0 = 10.0;     /* Starting energy in MeV (center-of-mass) */
    double estep = 10.0;  /* Energy step in MeV */

    /* Woods-Saxon potential parameters */
    double an = 0.5803;   /* Diffuseness in fm */
    double rn = 1.1132 * (pow(a1, 1.0/3) + pow(a2, 1.0/3));  /* Radius in fm */

    printf("\nParameters:\n");
    printf("  Angular momentum L = %d\n", l);
    printf("  Lagrange mesh: nr=%d, ns=%d, total points=%d\n", nr, ns, nr*ns);
    printf("  Channel radius: rmax = %.1f fm\n", rmax);

    /* ========================================================================
     * Step 1: Initialize R-matrix calculation
     *
     * hprmat_init(nr, ns, rmax, zrma) computes:
     *   - Lagrange-Legendre mesh points
     *   - Kinetic energy matrix elements
     *   - Bloch operator matrix elements
     *
     * Output: zrma = array of radial mesh points (size: nr*ns)
     * ======================================================================== */
    int ntot = nr * ns;
    double *zrma = (double*)malloc(ntot * sizeof(double));
    hprmat_init(nr, ns, rmax, zrma);

    printf("  Mesh range: %.4f to %.4f fm\n", zrma[0], zrma[ntot-1]);

    /* ========================================================================
     * Step 2: Build the coupling potential matrix
     *
     * cpot[ir + i*ntot + j*ntot*nch] = V_ij(r_ir) / (hbar^2/2mu)
     *
     * For single channel (nch=1): only cpot[ir] is used (column-major order)
     * The potential must be divided by hbar^2/(2*mu) for HPRMAT convention
     *
     * Here we use complex Woods-Saxon + Coulomb:
     *   V(r) = V_WS(r) + V_C(r)
     *   V_WS = -(V0 + i*W0) / (1 + exp((r-R)/a))  with V0=100, W0=10 MeV
     *   V_C  = Z1*Z2*e^2 / r
     * ======================================================================== */
    int nch = 1;  /* Number of channels (single channel for elastic scattering) */
    double _Complex *cpot = (double _Complex*)calloc(ntot * nch * nch,
                                                      sizeof(double _Complex));

    for (int i = 0; i < ntot; i++) {
        double r = zrma[i];
        /* Complex Woods-Saxon potential (absorptive) */
        double xx = 1.0 + exp((r - rn) / an);
        double _Complex cvn = (-100.0 - 10.0*I) / xx;  /* V0=100, W0=10 MeV */
        /* Point Coulomb potential */
        double vc = ze / r;
        /* Store potential divided by hbar^2/(2*mu) */
        /* For column-major order: cpot[ir, 0, 0] = cpot[ir] */
        cpot[i] = (cvn + vc) / hm;
    }

    /* ========================================================================
     * Step 3: Set up channel quantum numbers
     *
     * lval[i] = orbital angular momentum for channel i
     * For single channel elastic scattering, lval = {L}
     * ======================================================================== */
    int lval[1] = {l};

    printf("\nEnergy scan: %d points from %.1f to %.1f MeV\n",
           ne, e0, e0 + (ne-1)*estep);
    printf("------------------------------------------------------------\n");
    printf("%10s %14s %14s %10s\n", "E (MeV)", "Re(S)", "Im(S)", "|S|");
    printf("------------------------------------------------------------\n");

    /* Output arrays */
    double _Complex *cu = (double _Complex*)calloc(nch * nch,
                                                    sizeof(double _Complex));
    int nopen;

    /* ========================================================================
     * Step 4: Loop over energies and solver types
     *
     * For each energy E:
     *   - Calculate wave number: k = sqrt(2*mu*E) / hbar = sqrt(E/hm)
     *   - Calculate Sommerfeld parameter: eta = Z1*Z2*e^2*mu / (hbar^2*k)
     *   - Call hprmat_solve() to get the collision (S) matrix
     *
     * Solver types (defined in hprmat.h):
     *   HPRMAT_SOLVER_DENSE (1)    - Dense LAPACK, reference accuracy
     *   HPRMAT_SOLVER_MIXED (2)    - Mixed precision, faster
     *   HPRMAT_SOLVER_WOODBURY (3) - CPU optimized
     *   HPRMAT_SOLVER_GPU (4)      - GPU accelerated (if available)
     * ======================================================================== */
    int solvers[] = {HPRMAT_SOLVER_DENSE, HPRMAT_SOLVER_MIXED, HPRMAT_SOLVER_WOODBURY};
    const char *solver_names[] = {"Dense LAPACK", "Mixed Precision", "Woodbury"};

    for (int is = 0; is < 3; is++) {
        int solver = solvers[is];
        printf("\n--- Solver %d: %s ---\n", solver, solver_names[is]);

        for (int ie = 0; ie < ne; ie++) {
            double ecm = e0 + ie * estep;  /* Center-of-mass energy in MeV */

            /* Wave number k in fm^-1 */
            double qk_val = sqrt(ecm / hm);
            /* Sommerfeld parameter (dimensionless) */
            double eta_val = eta0 / qk_val;

            /* Arrays for HPRMAT (size = number of channels) */
            double qk[1] = {qk_val};
            double eta[1] = {eta_val};

            /* =================================================================
             * Call the R-matrix solver
             *
             * hprmat_solve(nch, lval, qk, eta, rmax, nr, ns,
             *              cpot, cpot_dim1, cu, &nopen, solver)
             *
             * Arguments:
             *   nch       - number of channels
             *   lval      - angular momentum array
             *   qk        - wave number array (positive=open, negative=closed)
             *   eta       - Sommerfeld parameter array
             *   rmax      - channel radius
             *   nr, ns    - mesh parameters
             *   cpot      - coupling potential matrix (column-major)
             *   cpot_dim1 - first dimension of cpot (= ntot)
             *   cu        - output collision (S) matrix
             *   nopen     - output: number of open channels
             *   solver    - solver type (1, 2, 3, or 4)
             * ================================================================= */
            hprmat_solve(nch, lval, qk, eta, rmax, nr, ns,
                         cpot, ntot, cu, &nopen, solver);

            /* =================================================================
             * Extract results
             *
             * For elastic scattering, S-matrix element S_11 gives:
             *   |S| = 1 for no absorption
             *   |S| < 1 indicates absorption (flux loss)
             *   Phase shift: delta = carg(S) / 2
             * ================================================================= */
            double _Complex S11 = cu[0];  /* cu[0,0] in column-major */
            printf("%10.3f %14.6e %14.6e %10.6f\n",
                   ecm, creal(S11), cimag(S11), cabs(S11));
        }
    }

    printf("------------------------------------------------------------\n");
    printf("\nTest completed!\n");

    /* Clean up */
    free(zrma);
    free(cpot);
    free(cu);

    return 0;
}
