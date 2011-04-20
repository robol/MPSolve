/***********************************************************
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.2, May 2001                  **
 **                                                        **
 **                      Written by                        **
 **       Dario Andrea Bini and Giuseppe Fiorentino        **
 **       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
 **                                                        **
 ** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
 ***********************************************************/

#include <mps/mps.h>

/**
 * @brief Compute Aberth correction for j-th root, without
 * selective correction.
 */
void
mps_faberth(mps_status* s, int j, cplx_t abcorr) {
    int i;
    cplx_t z;

    cplx_set(abcorr, cplx_zero);
    for (i = 0; i < s->n; i++) {
        if (i == j)
            continue;
        cplx_sub(z, s->froot[j], s->froot[i]);
        cplx_inv_eq(z);
        cplx_add_eq(abcorr, z);
    }
}

/**
 * @brief Compute Aberth correction for j-th root, without
 * selective correction.
 */
void
mps_daberth(mps_status* s, int j, cdpe_t abcorr) {
    int i;
    cdpe_t z;

    cdpe_set(abcorr, cdpe_zero);
    for (i = 0; i < s->n; i++) {
        if (i == j)
            continue;
        cdpe_sub(z, s->droot[j], s->droot[i]);
        cdpe_inv_eq(z);
        cdpe_add_eq(abcorr, z);
    }
}

/**
 * @brief Compute Aberth correction for j-th root, without
 * selective correction.
 */
void
mps_maberth(mps_status* s, int j, mpc_t abcorr) {
    int i;
    cdpe_t z, temp;
    tmpc_t diff;

    tmpc_init2(diff, s->mpwp);

    cdpe_set(temp, cdpe_zero);
    for (i = 0; i < s->n; i++) {
        if (i == j)
            continue;
        mpc_sub(diff, s->mroot[j], s->mroot[i]);
        mpc_get_cdpe(z, diff);
        cdpe_inv_eq(z);
        cdpe_add_eq(temp, z);
    }
    mpc_set_cdpe(abcorr, temp);

    tmpc_clear(diff);
}

/**
 * @brief Compute Aberth correction for the j-th root,
 * but only with other roots of the <code>jc</code>-th
 * cluster.
 */
void
mps_faberth_s(mps_status* s, int j, int jc, cplx_t abcorr) {
    int i, k;
    cplx_t z;

    cplx_set(abcorr, cplx_zero);
    for (i = s->punt[jc]; i < s->punt[jc + 1]; i++) {
        k = s->clust[i];
        if (k == j)
            continue;
        cplx_sub(z, s->froot[j], s->froot[k]);
        cplx_inv_eq(z);
        cplx_add_eq(abcorr, z);
    }
}

/**
 * @brief Compute Aberth correction for the j-th root,
 * but only with other roots of the <code>jc</code>-th
 * cluster.
 */
void
mps_daberth_s(mps_status* s, int j, int jc, cdpe_t abcorr) {
    int i, k;
    cdpe_t z;

    cdpe_set(abcorr, cdpe_zero);
    for (i = s->punt[jc]; i < s->punt[jc + 1]; i++) {
        k = s->clust[i];
        if (k == j)
            continue;
        cdpe_sub(z, s->droot[j], s->droot[k]);
        cdpe_inv_eq(z);
        cdpe_add_eq(abcorr, z);
    }
}

/**
 * @brief Compute Aberth correction for the j-th root,
 * but only with other roots of the <code>jc</code>-th
 * cluster.
 */
void
mps_maberth_s(mps_status* s, int j, int jc, mpc_t abcorr) {
    int i, k;
    cdpe_t z, temp;
    tmpc_t diff;

    tmpc_init2(diff, s->mpwp);

    cdpe_set(temp, cdpe_zero);
    for (i = s->punt[jc]; i < s->punt[jc + 1]; i++) {
        k = s->clust[i];
        if (k == j)
            continue;
        mpc_sub(diff, s->mroot[j], s->mroot[k]);
        mpc_get_cdpe(z, diff);
        cdpe_inv_eq(z);
        cdpe_add_eq(temp, z);
    }
    mpc_set_cdpe(abcorr, temp);

    tmpc_clear(diff);
}
