/*
 * rpca_api.h - thin, dependency-free C entry points for the Python binding.
 *
 * These wrap the proven RPCA numerical core (src/linear_algebra.c) behind a
 * buffer-in / buffer-out interface that is trivial to call from Cython/cffi/
 * ctypes: every output is written into a caller-allocated array, so the
 * binding never has to free C-owned memory.
 *
 * No GROMACS, no MDAnalysis - only LAPACK/BLAS (via the core) and libm.
 */
#ifndef RPCA_API_H
#define RPCA_API_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Simultaneous diagonalization of two covariance matrices followed by the
 * relative-PCA Kullback-Leibler ranking (Ahmad et al., JCTC 2019).
 *
 * Finds the generalized eigenbasis G such that G^T A G = I and
 * G^T B G = diag(geigval), then scores each mode by the KL divergence of
 * state B from state A and returns the modes ordered by descending KL.
 *
 * Inputs:
 *   n        - dimension (= 3 * n_atoms)
 *   A, B     - n*n symmetric covariance matrices (state A / state B).
 *              Row- or column-major is equivalent because they are symmetric.
 *              Neither is modified.
 *   mean_a,  - length-n mean coordinate vectors. Either may be NULL (treated
 *   mean_b     as zero), in which case the mean-shift term vanishes.
 *   algo     - 0: standard, 1: mean-fluctuation subspacing.
 *   verbose  - passed through to the core (0 = silent).
 *
 * Outputs (caller-allocated; size for the worst case r == n):
 *   geigval  - length >= n; first r entries are the generalized eigenvalues.
 *   gevec    - length >= n*n; first n*r entries are the eigenvectors in
 *              COLUMN-MAJOR n x r layout (mode k = gevec[k*n .. k*n+n-1]).
 *   kl       - length >= n; per-mode KL divergence, descending.
 *   kl_m     - length >= n; per-mode KL contribution from the mean shift.
 *   acc_kl   - length >= n; cumulative percentage of total KL.
 *   rank_out - the effective rank r (<= n). Only the first r entries / columns
 *              of the output buffers are meaningful.
 *
 * Returns 0 on success, non-zero on error.
 */
int rpca_sdiag(int n,
               const double *A, const double *B,
               const double *mean_a, const double *mean_b,
               int algo, int verbose,
               double *geigval, double *gevec,
               double *kl, double *kl_m, double *acc_kl,
               int *rank_out);

#ifdef __cplusplus
}
#endif

#endif /* RPCA_API_H */
