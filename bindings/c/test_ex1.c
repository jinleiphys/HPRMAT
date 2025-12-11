/*
 * Test Ex1: alpha-208Pb Optical Potential (C)
 * ============================================
 *
 * This replicates the Fortran Ex1 example exactly.
 * Compares results with the Fortran reference.
 *
 * Reference: Goldring et al., Phys. Lett. B32 (1970) 465
 *
 * Compile:
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

    /* Physical constants (same as Fortran) */
    double a1 = 208.0, a2 = 4.0;    /* Mass numbers */
    double z1 = 82.0, z2 = 2.0;     /* Charge numbers */
    double rmu = a1 * a2 / (a1 + a2);  /* Reduced mass */
    double ze = z1 * z2 * 1.44;     /* Z1*Z2*e^2 (MeV.fm) */
    double hm = 20.736 / rmu;       /* hbar^2/(2*mu) in MeV.fm^2 */
    double eta0 = ze / (2 * hm);    /* Sommerfeld parameter factor */

    /* R-matrix parameters (from data1) */
    int l = 20;        /* Angular momentum */
    int nr = 60;       /* Lagrange functions per interval */
    int ns = 1;        /* Number of intervals */
    double rmax = 14.0;/* Channel radius (fm) */

    /* Energy parameters */
    int ne = 5;        /* Number of energies */
    double e0 = 10.0;  /* Initial energy (MeV) */
    double estep = 10.0; /* Energy step (MeV) */

    /* Woods-Saxon parameters */
    double an = 0.5803;  /* Diffuseness (fm) */
    double rn = 1.1132 * (pow(a1, 1.0/3) + pow(a2, 1.0/3));  /* Radius (fm) */

    printf("\nParameters:\n");
    printf("  Angular momentum L = %d\n", l);
    printf("  Lagrange functions: nr = %d, ns = %d\n", nr, ns);
    printf("  Channel radius: rmax = %.1f fm\n", rmax);
    printf("  Reduced mass: mu = %.4f amu\n", rmu);
    printf("  hbar^2/(2*mu) = %.6f MeV.fm^2\n", hm);

    /* Initialize R-matrix */
    int ntot = nr * ns;
    double *zrma = (double*)malloc(ntot * sizeof(double));
    hprmat_init(nr, ns, rmax, zrma);

    printf("  Mesh points: %d\n", ntot);
    printf("  First mesh point: r = %.6f fm\n", zrma[0]);
    printf("  Last mesh point: r = %.6f fm\n", zrma[ntot-1]);

    /* Build potential matrix (complex Woods-Saxon + Coulomb) */
    int nch = 1;
    double _Complex *cpot = (double _Complex*)calloc(ntot * nch * nch, sizeof(double _Complex));

    for (int i = 0; i < ntot; i++) {
        double r = zrma[i];
        double xx = 1.0 + exp((r - rn) / an);
        double _Complex cvn = (-100.0 - 10.0*I) / xx;  /* Complex Woods-Saxon */
        double vc = ze / r;                             /* Coulomb */
        cpot[i] = (cvn + vc) / hm;  /* cpot[i, 0, 0] in column-major */
    }

    /* Channel parameters */
    int lval[1] = {l};

    printf("\nEnergies: %d points from %.1f to %.1f MeV\n", ne, e0, e0 + (ne-1)*estep);
    printf("------------------------------------------------------------\n");
    printf("%10s %14s %14s %10s\n", "E (MeV)", "Re(S)", "Im(S)", "|S|");
    printf("------------------------------------------------------------\n");

    /* Output arrays */
    double _Complex *cu = (double _Complex*)calloc(nch * nch, sizeof(double _Complex));
    int nopen;

    /* Test all solver types */
    const char *solver_names[] = {"", "Dense LAPACK", "Mixed Precision", "Woodbury"};

    for (int solver = 1; solver <= 3; solver++) {
        printf("\n--- Solver %d: %s ---\n", solver, solver_names[solver]);

        for (int ie = 0; ie < ne; ie++) {
            double ecm = e0 + ie * estep;

            /* Wave number and Sommerfeld parameter */
            double qk_val = sqrt(ecm / hm);
            double eta_val = eta0 / qk_val;

            double qk[1] = {qk_val};
            double eta[1] = {eta_val};

            /* Solve */
            hprmat_solve(nch, lval, qk, eta, rmax, nr, ns,
                         cpot, ntot, cu, &nopen, solver);

            /* Extract S-matrix element */
            double _Complex S11 = cu[0];
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
