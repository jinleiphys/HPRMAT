/*
 * Alpha-Alpha Scattering Example (C)
 * ===================================
 *
 * This example demonstrates using the HPRMAT C interface for single-channel
 * alpha-alpha elastic scattering with the Ali-Bodmer potential.
 *
 * Compile with:
 *   gcc -o example_alpha_alpha example_alpha_alpha.c -L../../lib -lhprmat -lm
 *
 * Run with:
 *   LD_LIBRARY_PATH=../../lib ./example_alpha_alpha
 *   (or on macOS: DYLD_LIBRARY_PATH=../../lib ./example_alpha_alpha)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "hprmat.h"

/*
 * Ali-Bodmer alpha-alpha potential (deep version).
 * V(r) = 500 * exp(-0.7^2 * r^2) - 130 * exp(-0.475^2 * r^2)  (MeV)
 */
double ali_bodmer_potential(double r) {
    return 500.0 * exp(-0.7*0.7 * r*r) - 130.0 * exp(-0.475*0.475 * r*r);
}

int main() {
    /* Physical parameters */
    double mu = 4.0 * 938.0 / 2.0;  /* Reduced mass for alpha-alpha (MeV/c^2) */
    double hbarc = 197.3;            /* MeV.fm */

    /* R-matrix parameters */
    int nr = 30;        /* Lagrange functions */
    int ns = 1;         /* Single interval */
    double rmax = 10.0; /* Channel radius (fm) */

    /* Channel parameters */
    int nch = 1;
    int lval[1] = {0};      /* s-wave */
    double eta[1] = {0.0};  /* No Coulomb */

    /* Allocate mesh array */
    double *zrma = (double*)malloc(nr * ns * sizeof(double));

    /* Initialize R-matrix solver */
    printf("Initializing R-matrix solver...\n");
    hprmat_init(nr, ns, rmax, zrma);
    hprmat_set_solver(HPRMAT_SOLVER_DENSE);

    char version[64];
    hprmat_version(version, 64);
    printf("  Version: %s\n", version);
    printf("  Mesh points: %d\n", nr * ns);

    /* Build potential matrix on Lagrange mesh */
    double _Complex *cpot = (double _Complex*)calloc(nr * ns * nch * nch,
                                                      sizeof(double _Complex));

    for (int ir = 0; ir < nr * ns; ir++) {
        double r = zrma[ir];
        double V = 2.0 * mu / (hbarc * hbarc) * ali_bodmer_potential(r);
        cpot[ir] = V + 0.0*I;  /* cpot[ir, 0, 0] in column-major */
    }

    /* Output collision matrix */
    double _Complex *cu = (double _Complex*)calloc(nch * nch, sizeof(double _Complex));
    int nopen;

    /* Energy scan */
    double energies[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    int n_energies = 5;

    printf("\nAlpha-Alpha Elastic Scattering (Ali-Bodmer Potential)\n");
    printf("============================================================\n");
    printf("%12s %12s %12s %12s %12s\n", "E_cm (MeV)", "k (fm^-1)", "Re(S)", "Im(S)", "delta (deg)");
    printf("------------------------------------------------------------\n");

    for (int ie = 0; ie < n_energies; ie++) {
        double E = energies[ie];

        /* Wave number */
        double k = sqrt(2.0 * mu * E) / hbarc;
        double qk[1] = {k};

        /* Solve */
        hprmat_solve(nch, lval, qk, eta, rmax, nr, ns,
                     cpot, nr * ns, cu, &nopen, HPRMAT_SOLVER_DENSE);

        /* Extract phase shift */
        double _Complex S_11 = cu[0];  /* cu[0, 0] */
        double delta = carg(S_11) / 2.0 * 180.0 / M_PI;

        printf("%12.2f %12.4f %12.6f %12.6f %12.2f\n",
               E, k, creal(S_11), cimag(S_11), delta);
    }

    /* Compare solvers at E = 4 MeV */
    printf("\n\nSolver Comparison at E = 4 MeV\n");
    printf("============================================================\n");

    double E = 4.0;
    double k = sqrt(2.0 * mu * E) / hbarc;
    double qk[1] = {k};

    int solvers[] = {HPRMAT_SOLVER_DENSE, HPRMAT_SOLVER_MIXED, HPRMAT_SOLVER_WOODBURY};
    const char *solver_names[] = {"Dense LAPACK", "Mixed Precision", "Woodbury"};
    int n_solvers = 3;

    double _Complex results[3];

    for (int is = 0; is < n_solvers; is++) {
        char info[128];
        hprmat_solver_info(solvers[is], info, 128);

        hprmat_solve(nch, lval, qk, eta, rmax, nr, ns,
                     cpot, nr * ns, cu, &nopen, solvers[is]);

        results[is] = cu[0];
        printf("%20s: S = %12.8f + %12.8fi\n",
               solver_names[is], creal(cu[0]), cimag(cu[0]));
    }

    /* Show differences */
    printf("\nDifference from Dense LAPACK:\n");
    for (int is = 1; is < n_solvers; is++) {
        double diff = cabs(results[is] - results[0]);
        printf("  %20s: |dS| = %.2e\n", solver_names[is], diff);
    }

    /* Clean up */
    free(zrma);
    free(cpot);
    free(cu);

    printf("\nDone!\n");
    return 0;
}
