# cython: language_level=3
"""Cython binding for the RPCA simultaneous-diagonalization core."""

import numpy as np
cimport numpy as cnp

cdef extern from "rpca_api.h":
    int rpca_sdiag(int n,
                   const double *A, const double *B,
                   const double *mean_a, const double *mean_b,
                   int algo, int verbose,
                   double *geigval, double *gevec,
                   double *kl, double *kl_m, double *acc_kl,
                   int *rank_out) nogil

cnp.import_array()


def sdiag(cov_a, cov_b, mean_a=None, mean_b=None, int algo=0, int verbose=0):
    """Simultaneous diagonalization + relative-PCA KL ranking of two states.

    Parameters
    ----------
    cov_a, cov_b : (n, n) array_like
        Symmetric covariance matrices of state A and state B.
    mean_a, mean_b : (n,) array_like, optional
        Mean coordinate vectors. If omitted the mean-shift term is zero.
    algo : int
        0 = standard diagonalization, 1 = mean-fluctuation subspacing.
    verbose : int
        Verbosity passed to the C core.

    Returns
    -------
    dict with keys ``geigval``, ``gevec`` (n x rank, column k = mode k),
    ``kl``, ``kl_m``, ``acc_kl`` (all length ``rank``), ``rank``,
    ``sum_kl`` and ``sum_kl_m``.
    """
    cdef cnp.ndarray[double, ndim=2, mode="c"] A = np.ascontiguousarray(cov_a, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2, mode="c"] B = np.ascontiguousarray(cov_b, dtype=np.float64)
    cdef int n = A.shape[0]

    if A.shape[1] != n or B.shape[0] != n or B.shape[1] != n:
        raise ValueError("cov_a and cov_b must be square matrices of the same size")

    cdef cnp.ndarray[double, ndim=1, mode="c"] ma
    cdef cnp.ndarray[double, ndim=1, mode="c"] mb
    ma = (np.zeros(n, dtype=np.float64) if mean_a is None
          else np.ascontiguousarray(mean_a, dtype=np.float64).ravel())
    mb = (np.zeros(n, dtype=np.float64) if mean_b is None
          else np.ascontiguousarray(mean_b, dtype=np.float64).ravel())
    if ma.shape[0] != n or mb.shape[0] != n:
        raise ValueError("mean vectors must have length n")

    cdef cnp.ndarray[double, ndim=1, mode="c"] geigval = np.zeros(n, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=1, mode="c"] gevec   = np.zeros(n * n, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=1, mode="c"] kl      = np.zeros(n, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=1, mode="c"] kl_m    = np.zeros(n, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=1, mode="c"] acc_kl  = np.zeros(n, dtype=np.float64)
    cdef int rank = 0
    cdef int status

    with nogil:
        status = rpca_sdiag(n, &A[0, 0], &B[0, 0], &ma[0], &mb[0],
                            algo, verbose,
                            &geigval[0], &gevec[0],
                            &kl[0], &kl_m[0], &acc_kl[0], &rank)
    if status != 0:
        raise RuntimeError(f"rpca_sdiag failed with status {status}")

    r = rank
    # Core writes column-major n x r: column k occupies gevec[k*n : (k+1)*n].
    gevec_2d = np.ascontiguousarray(gevec[:n * r].reshape((r, n)).T)

    return {
        "geigval": geigval[:r].copy(),
        "gevec": gevec_2d,
        "kl": kl[:r].copy(),
        "kl_m": kl_m[:r].copy(),
        "acc_kl": acc_kl[:r].copy(),
        "rank": r,
        "sum_kl": float(kl[:r].sum()),
        "sum_kl_m": float(kl_m[:r].sum()),
    }
