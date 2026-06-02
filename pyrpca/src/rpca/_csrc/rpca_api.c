/*
 * rpca_api.c - implementation of the thin C entry points.
 *
 * The heavy lifting (whitening transformation + LAPACK eigendecomposition) is
 * done by Simultaneous_Diagonalization() in the existing core; this file only
 * adds the lightweight KL scoring and ordering glue and marshals the result
 * into caller-owned buffers.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "linear_algebra.h"
#include "rpca_api.h"

/* Indices that sort `key` in descending order. r is the rank of a covariance
 * matrix (small relative to a trajectory), so an insertion sort is plenty. */
static void argsort_desc(const double *key, int r, int *idx)
{
    int i, j, tmp;
    for (i = 0; i < r; i++) idx[i] = i;
    for (i = 1; i < r; i++) {
        tmp = idx[i];
        j = i - 1;
        while (j >= 0 && key[idx[j]] < key[tmp]) {
            idx[j + 1] = idx[j];
            j--;
        }
        idx[j + 1] = tmp;
    }
}

int rpca_sdiag(int n,
               const double *A, const double *B,
               const double *mean_a, const double *mean_b,
               int algo, int verbose,
               double *geigval, double *gevec,
               double *kl, double *kl_m, double *acc_kl,
               int *rank_out)
{
    double *gv = NULL;   /* core-allocated eigenvalues, length r            */
    double *ge = NULL;   /* core-allocated eigenvectors, column-major n x r */
    double *dav = NULL;  /* mean shift (state B relative to state A)         */
    double *kl_tmp = NULL, *klm_tmp = NULL;
    int *order = NULL;
    int r = 0, i, k, rc = 0;

    if (n <= 0 || A == NULL || B == NULL || rank_out == NULL) return 1;

    dav = (double *) malloc((size_t) n * sizeof(double));
    if (!dav) { rc = 2; goto cleanup; }
    for (i = 0; i < n; i++) {
        double ma = mean_a ? mean_a[i] : 0.0;
        double mb = mean_b ? mean_b[i] : 0.0;
        dav[i] = mb - ma;
    }

    /* The core takes non-const pointers but copies A/B internally before any
     * destructive LAPACK call, so the caller's matrices are left untouched. */
    if (algo == 1) {
        Simultaneous_Diagonalization_subspacing(n, (double *) A, (double *) B,
                                                &gv, &ge, &r, dav, verbose);
    } else {
        Simultaneous_Diagonalization(n, (double *) A, (double *) B,
                                     &gv, &ge, &r, verbose);
    }
    if (r <= 0 || gv == NULL || ge == NULL) { rc = 3; goto cleanup; }

    kl_tmp  = (double *) malloc((size_t) r * sizeof(double));
    klm_tmp = (double *) malloc((size_t) r * sizeof(double));
    order   = (int *)    malloc((size_t) r * sizeof(int));
    if (!kl_tmp || !klm_tmp || !order) { rc = 4; goto cleanup; }

    /* Per-mode KL divergence:
     *   kl_m = 0.5 * (g_k . (mean_b - mean_a))^2
     *   kl   = 0.5 * (geigval_k - log(geigval_k) - 1) + kl_m            */
    for (k = 0; k < r; k++) {
        const double *vk = ge + (size_t) k * n;   /* column k */
        double proj = 0.0;
        for (i = 0; i < n; i++) proj += vk[i] * dav[i];
        klm_tmp[k] = 0.5 * proj * proj;
        kl_tmp[k]  = 0.5 * (gv[k] - log(gv[k]) - 1.0) + klm_tmp[k];
    }

    argsort_desc(kl_tmp, r, order);

    double kl_sum = 0.0;
    for (k = 0; k < r; k++) kl_sum += kl_tmp[k];

    double run = 0.0;
    for (k = 0; k < r; k++) {
        int src = order[k];
        geigval[k] = gv[src];
        kl[k]      = kl_tmp[src];
        kl_m[k]    = klm_tmp[src];
        run += kl_tmp[src];
        acc_kl[k]  = (kl_sum != 0.0) ? (run / kl_sum * 100.0) : 0.0;
        memcpy(gevec + (size_t) k * n, ge + (size_t) src * n,
               (size_t) n * sizeof(double));
    }

    *rank_out = r;

cleanup:
    free(dav);
    free(kl_tmp);
    free(klm_tmp);
    free(order);
    free(gv);   /* allocated by the core with malloc */
    free(ge);
    return rc;
}
