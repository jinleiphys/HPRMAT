"""
HPRMAT - High-Performance R-Matrix Solver
==========================================

Python interface for the HPRMAT Fortran library using f2py.

Build the Fortran extension first:
    cd bindings/python
    make

Example
-------
>>> import hprmat
>>> import numpy as np
>>>
>>> # Initialize
>>> nr, ns, rmax = 20, 1, 10.0
>>> solver = hprmat.RMatrixSolver(nr, ns, rmax)
>>>
>>> # Set up a single channel problem
>>> nch = 1
>>> lval = np.array([0], dtype=np.int32)
>>> qk = np.array([0.5], dtype=np.float64)
>>> eta = np.array([0.0], dtype=np.float64)
>>>
>>> # Build potential
>>> cpot = np.zeros((nr*ns, nch, nch), dtype=np.complex128, order='F')
>>> for ir, r in enumerate(solver.mesh):
...     cpot[ir, 0, 0] = -50.0 * np.exp(-r**2 / 4.0)
>>>
>>> # Solve
>>> S_matrix, nopen = solver.solve(lval, qk, eta, cpot)
>>> print(f"S-matrix: {S_matrix[0,0]:.6f}")

Author: Jin Lei
Date: December 2025
"""

import numpy as np

# Try to import the f2py-generated module
try:
    from . import hprmat_fortran as _f
    _HAS_FORTRAN = True
except ImportError:
    try:
        import hprmat_fortran as _f
        _HAS_FORTRAN = True
    except ImportError:
        _HAS_FORTRAN = False
        _f = None

__version__ = "1.0.0"
__author__ = "Jin Lei"

# Solver type constants
SOLVER_DENSE = 1      # Dense LAPACK ZGESV (reference)
SOLVER_MIXED = 2      # Mixed precision
SOLVER_WOODBURY = 3   # Woodbury-Kinetic (CPU optimized)
SOLVER_GPU = 4        # GPU cuSOLVER

SOLVER_NAMES = {
    SOLVER_DENSE: "Dense LAPACK ZGESV (reference)",
    SOLVER_MIXED: "Mixed Precision (single + double refinement)",
    SOLVER_WOODBURY: "Woodbury-Kinetic (CPU optimized)",
    SOLVER_GPU: "GPU cuSOLVER (NVIDIA GPU)",
}


def check_fortran():
    """Check if Fortran extension is available."""
    if not _HAS_FORTRAN:
        raise ImportError(
            "HPRMAT Fortran extension not found. Build it first:\n"
            "  cd bindings/python && make"
        )


class RMatrixSolver:
    """
    High-performance R-matrix solver for nuclear scattering.

    Parameters
    ----------
    nr : int
        Number of Lagrange functions per interval (typically 20-40)
    ns : int
        Number of intervals (typically 1)
    rmax : float
        R-matrix channel radius in fm
    solver : int, optional
        Solver type (default: SOLVER_DENSE=1)
        - 1: Dense LAPACK, most accurate
        - 2: Mixed precision, faster
        - 3: CPU optimized (Woodbury)
        - 4: GPU accelerated (requires CUDA)

    Attributes
    ----------
    mesh : ndarray
        Lagrange mesh abscissas (size: nr*ns)
    nr, ns, rmax : int, int, float
        Problem parameters

    Examples
    --------
    >>> solver = RMatrixSolver(nr=30, ns=1, rmax=12.0)
    >>> print(f"Mesh points: {len(solver.mesh)}")
    """

    def __init__(self, nr, ns, rmax, solver=SOLVER_DENSE):
        check_fortran()

        self.nr = int(nr)
        self.ns = int(ns)
        self.rmax = float(rmax)
        self.solver = int(solver)

        # Initialize and get mesh
        self.mesh = _f.py_rmat_init(self.nr, self.ns, self.rmax)

        # Set solver type
        _f.py_set_solver(self.solver)

    def solve(self, lval, qk, eta, cpot, solver=None):
        """
        Solve the R-matrix scattering problem.

        Parameters
        ----------
        lval : array_like
            Angular momentum for each channel (size: nch)
        qk : array_like
            Wave numbers (size: nch), positive=open, negative=closed
        eta : array_like
            Sommerfeld parameters (size: nch)
        cpot : array_like
            Coupling potential (size: nr*ns, nch, nch), Fortran order
        solver : int, optional
            Override solver type for this calculation

        Returns
        -------
        cu : ndarray
            Collision (S) matrix (size: nch, nch)
        nopen : int
            Number of open channels
        """
        # Ensure correct types and memory layout
        lval = np.asarray(lval, dtype=np.int32)
        qk = np.asarray(qk, dtype=np.float64)
        eta = np.asarray(eta, dtype=np.float64)
        cpot = np.asfortranarray(cpot, dtype=np.complex128)

        nch = len(lval)
        use_solver = solver if solver is not None else self.solver

        cu, nopen = _f.py_rmatrix(
            nch, lval, qk, eta, self.rmax, self.nr, self.ns,
            cpot, use_solver
        )

        return cu, nopen

    def solve_with_wavefunction(self, lval, qk, eta, cpot,
                                 entrance_channels=None, solver=None):
        """
        Solve R-matrix problem and compute wave functions.

        Parameters
        ----------
        lval, qk, eta, cpot : array_like
            Same as solve()
        entrance_channels : array_like, optional
            Indices of entrance channels (1-based). Default: [1]
        solver : int, optional
            Override solver type

        Returns
        -------
        cu : ndarray
            Collision matrix (nch, nch)
        cf : ndarray
            Wave functions (nr*ns, nch, n_entrance)
        nopen : int
            Number of open channels
        """
        lval = np.asarray(lval, dtype=np.int32)
        qk = np.asarray(qk, dtype=np.float64)
        eta = np.asarray(eta, dtype=np.float64)
        cpot = np.asfortranarray(cpot, dtype=np.complex128)

        nch = len(lval)
        use_solver = solver if solver is not None else self.solver

        if entrance_channels is None:
            nvc = np.array([1], dtype=np.int32)
        else:
            nvc = np.asarray(entrance_channels, dtype=np.int32)

        nc_entrance = len(nvc)

        cu, cf, nopen = _f.py_rmatrix_wf(
            nch, lval, qk, eta, self.rmax, self.nr, self.ns,
            cpot, nc_entrance, nvc, use_solver
        )

        return cu, cf, nopen

    def solve_nonlocal(self, lval, qk, eta, cpot, cpnl,
                        entrance_channels=None, solver=None):
        """
        Solve R-matrix problem with non-local potential.

        Parameters
        ----------
        lval, qk, eta, cpot : array_like
            Same as solve()
        cpnl : array_like
            Non-local potential (size: nr*nr, nch, nch)
        entrance_channels : array_like, optional
            Entrance channels (1-based)
        solver : int, optional
            Override solver type

        Returns
        -------
        cu, cf, nopen : ndarray, ndarray, int
            Same as solve_with_wavefunction()
        """
        lval = np.asarray(lval, dtype=np.int32)
        qk = np.asarray(qk, dtype=np.float64)
        eta = np.asarray(eta, dtype=np.float64)
        cpot = np.asfortranarray(cpot, dtype=np.complex128)
        cpnl = np.asfortranarray(cpnl, dtype=np.complex128)

        nch = len(lval)
        use_solver = solver if solver is not None else self.solver

        if entrance_channels is None:
            nvc = np.array([1], dtype=np.int32)
        else:
            nvc = np.asarray(entrance_channels, dtype=np.int32)

        nc_entrance = len(nvc)

        cu, cf, nopen = _f.py_rmatrix_nonlocal(
            nch, lval, qk, eta, self.rmax, self.nr, self.ns,
            cpot, cpnl, nc_entrance, nvc, use_solver
        )

        return cu, cf, nopen

    def interpolate_wavefunction(self, lval, qk, eta, cu, cf, nopen,
                                  channel, entrance, npoin, h=None):
        """
        Interpolate wave function onto uniform mesh.

        Parameters
        ----------
        lval, qk, eta : array_like
            Channel parameters
        cu : ndarray
            Collision matrix from solve_with_wavefunction
        cf : ndarray
            Wave functions from solve_with_wavefunction
        nopen : int
            Number of open channels
        channel : int
            Channel index (1-based)
        entrance : int
            Entrance channel number (1-based)
        npoin : int
            Number of output points
        h : float, optional
            Step size (default: 2*rmax/npoin)

        Returns
        -------
        r : ndarray
            Radial coordinates
        wf : ndarray
            Wave function values (complex)
        """
        lval = np.asarray(lval, dtype=np.int32)
        qk = np.asarray(qk, dtype=np.float64)
        eta = np.asarray(eta, dtype=np.float64)
        cu = np.asfortranarray(cu, dtype=np.complex128)
        cf = np.asfortranarray(cf, dtype=np.complex128)

        nch = len(lval)
        nom = cf.shape[2] if cf.ndim == 3 else 1

        if h is None:
            h = 2.0 * self.rmax / npoin

        cwftab = _f.py_wf_print(
            nch, lval, qk, eta, self.rmax, self.nr, self.ns,
            cu, nopen, cf, self.mesh, channel, nom, npoin, h
        )

        r = np.arange(1, npoin + 1) * h
        return r, cwftab

    def set_solver(self, solver):
        """Set the default solver type."""
        self.solver = int(solver)
        _f.py_set_solver(self.solver)

    @property
    def solver_name(self):
        """Get the name of the current solver."""
        return SOLVER_NAMES.get(self.solver, "Unknown")

    @staticmethod
    def available_solvers():
        """List available solvers."""
        return SOLVER_NAMES.copy()


# Convenience functions for direct solver access
def solve_dense(cmat, B_vector, nch, nlag, normfac=1.0):
    """Direct call to dense LAPACK solver."""
    check_fortran()
    cmat = np.asfortranarray(cmat, dtype=np.complex128)
    B_vector = np.asarray(B_vector, dtype=np.complex128)
    return _f.py_solve_rmatrix_dense(cmat, B_vector, nch, nlag, normfac)


def solve_mixed(cmat, B_vector, nch, nlag, normfac=1.0):
    """Direct call to mixed precision solver."""
    check_fortran()
    cmat = np.asfortranarray(cmat, dtype=np.complex128)
    B_vector = np.asarray(B_vector, dtype=np.complex128)
    return _f.py_solve_rmatrix_mixed(cmat, B_vector, nch, nlag, normfac)


def solve_woodbury(cmat, B_vector, nch, nlag, normfac=1.0):
    """Direct call to Woodbury solver."""
    check_fortran()
    cmat = np.asfortranarray(cmat, dtype=np.complex128)
    B_vector = np.asarray(B_vector, dtype=np.complex128)
    return _f.py_solve_rmatrix_woodbury(cmat, B_vector, nch, nlag, normfac)


def solve_gpu(cmat, B_vector, nch, nlag, normfac=1.0):
    """Direct call to GPU solver."""
    check_fortran()
    cmat = np.asfortranarray(cmat, dtype=np.complex128)
    B_vector = np.asarray(B_vector, dtype=np.complex128)
    return _f.py_solve_rmatrix_gpu(cmat, B_vector, nch, nlag, normfac)
