/*
 * HPRMAT - High-Performance R-Matrix Solver
 * C Interface Header
 *
 * This header provides C bindings for the HPRMAT Fortran library.
 * Compatible with C, C++, Rust, Go, and other languages.
 *
 * Author: Jin Lei
 * Date: December 2025
 */

#ifndef HPRMAT_H
#define HPRMAT_H

#include <complex.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Solver type constants */
#define HPRMAT_SOLVER_DENSE     1   /* Dense LAPACK ZGESV (reference) */
#define HPRMAT_SOLVER_MIXED     2   /* Mixed precision (single + double refinement) */
#define HPRMAT_SOLVER_WOODBURY  3   /* Woodbury-Kinetic (CPU optimized) */
#define HPRMAT_SOLVER_GPU       4   /* GPU cuSOLVER (NVIDIA GPU) */

/*
 * Initialize R-matrix calculation
 *
 * Must be called before hprmat_solve() to set up Lagrange mesh and
 * kinetic energy matrix elements.
 *
 * Parameters:
 *   nr   - Number of Lagrange functions per interval
 *   ns   - Number of intervals (typically 1)
 *   rmax - R-matrix channel radius in fm
 *   zrma - Output array for mesh abscissas (size: nr*ns)
 */
void hprmat_init(int nr, int ns, double rmax, double* zrma);

/*
 * Set the default solver type
 *
 * Parameters:
 *   solver - Solver type (1-4), see HPRMAT_SOLVER_* constants
 */
void hprmat_set_solver(int solver);

/*
 * Get the current default solver type
 *
 * Returns:
 *   Current solver type (1-4)
 */
int hprmat_get_solver(void);

/*
 * Main R-matrix solver (simplified interface)
 *
 * Computes the collision matrix for a given scattering problem.
 *
 * Parameters:
 *   nch       - Number of channels
 *   lval      - Angular momentum for each channel (size: nch)
 *   qk        - Wave numbers (size: nch), positive=open, negative=closed
 *   eta       - Sommerfeld parameters (size: nch)
 *   rmax      - Channel radius in fm
 *   nr        - Lagrange functions per interval
 *   ns        - Number of intervals
 *   cpot      - Coupling potentials, column-major (size: cpot_dim1 * nch * nch)
 *               cpot[ir + i*cpot_dim1 + j*cpot_dim1*nch] = V_ij(r_ir)
 *   cpot_dim1 - First dimension of cpot (typically nr*ns)
 *   cu        - Output collision matrix (size: nch * nch)
 *   nopen     - Output: number of open channels
 *   solver    - Solver type (1-4), or 0 for default
 */
void hprmat_solve(int nch, int* lval, double* qk, double* eta,
                  double rmax, int nr, int ns,
                  double _Complex* cpot, int cpot_dim1,
                  double _Complex* cu, int* nopen,
                  int solver);

/*
 * Full R-matrix solver with wave function output
 *
 * Parameters:
 *   nch         - Number of channels
 *   lval        - Angular momentum for each channel (size: nch)
 *   qk          - Wave numbers (size: nch)
 *   eta         - Sommerfeld parameters (size: nch)
 *   rmax        - Channel radius in fm
 *   nr          - Lagrange functions per interval
 *   ns          - Number of intervals
 *   cpot        - Coupling potentials (size: cpot_dim1 * nch * nch)
 *   cpot_dim1   - First dimension of cpot
 *   cu          - Output collision matrix (size: nch * nch)
 *   cf          - Output wave functions (size: cf_dim1 * cf_dim2 * nc_entrance)
 *   cf_dim1     - First dimension of cf (typically nr*ns)
 *   cf_dim2     - Second dimension of cf (typically nch)
 *   nc_entrance - Number of entrance channels
 *   nvc         - Entrance channel indices (size: nc_entrance)
 *   cpnl        - Non-local potentials (size: cpnl_dim1 * nch * nch), NULL if none
 *   cpnl_dim1   - First dimension of cpnl (0 if no non-local potential)
 *   nopen       - Output: number of open channels
 *   solver      - Solver type (1-4), or 0 for default
 */
void hprmat_solve_full(int nch, int* lval, double* qk, double* eta,
                       double rmax, int nr, int ns,
                       double _Complex* cpot, int cpot_dim1,
                       double _Complex* cu,
                       double _Complex* cf, int cf_dim1, int cf_dim2,
                       int nc_entrance, int* nvc,
                       double _Complex* cpnl, int cpnl_dim1,
                       int* nopen, int solver);

/*
 * Calculate wave function on uniform mesh
 *
 * Interpolates the internal wave function and calculates external asymptotics.
 *
 * Parameters:
 *   nch      - Number of channels
 *   lval     - Angular momentum values (size: nch)
 *   qk       - Wave numbers (size: nch)
 *   eta      - Sommerfeld parameters (size: nch)
 *   rmax     - Channel radius
 *   nr       - Lagrange functions per interval
 *   ns       - Number of intervals
 *   cu       - Collision matrix (size: nch * nch)
 *   nopen    - Number of open channels
 *   cf       - Wave functions from hprmat_solve_full (size: cf_dim1 * cf_dim2 * nom)
 *   cf_dim1  - First dimension of cf
 *   cf_dim2  - Second dimension of cf
 *   zrma     - Mesh abscissas from hprmat_init (size: nr*ns)
 *   iv       - Channel index for output
 *   nom      - Entrance channel number
 *   npoin    - Number of output points
 *   h        - Step size for output mesh
 *   cwftab   - Output wave function values (size: npoin)
 */
void hprmat_wavefunction(int nch, int* lval, double* qk, double* eta,
                         double rmax, int nr, int ns,
                         double _Complex* cu, int nopen,
                         double _Complex* cf, int cf_dim1, int cf_dim2,
                         double* zrma, int iv, int nom,
                         int npoin, double h,
                         double _Complex* cwftab);

/*
 * Low-level linear system solver
 *
 * Directly solves the R-matrix linear system (for advanced users).
 *
 * Parameters:
 *   cmat     - Coupling matrix (size: nch*nlag x nch*nlag)
 *   B_vector - Boundary vector (size: nlag)
 *   nch      - Number of channels
 *   nlag     - Number of Lagrange functions
 *   normfac  - Normalization factor
 *   Rmat     - Output R-matrix (size: nch x nch)
 *   solver   - Solver type (1-4)
 */
void hprmat_solve_linear(double _Complex* cmat, double _Complex* B_vector,
                         int nch, int nlag, double normfac,
                         double _Complex* Rmat, int solver);

/*
 * Get HPRMAT version string
 *
 * Parameters:
 *   version - Output buffer for version string
 *   maxlen  - Maximum length of version buffer
 */
void hprmat_version(char* version, int maxlen);

/*
 * Get solver description string
 *
 * Parameters:
 *   solver - Solver type (1-4)
 *   info   - Output buffer for description
 *   maxlen - Maximum length of info buffer
 */
void hprmat_solver_info(int solver, char* info, int maxlen);

#ifdef __cplusplus
}
#endif

#endif /* HPRMAT_H */
