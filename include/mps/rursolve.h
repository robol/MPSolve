/***********************************************************
**       Multiprecision Polynomial Solver (MPSolve)       **
**              Version 2.1, september 1999               **
**                                                        **
**                      Written by                        **
**       Dario Andrea Bini and Giuseppe Fiorentino        **
**       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
**                                                        **
** (C) 1999, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

#ifndef MPS_RURSOLVE_H
#define MPS_RURSOLVE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <mps/core.h>

extern mpz_t *mpdemo;	/* imaginary part of the integer input coeff. */

/* functions in main.c */
void mps_rursolve(mps_status* s);

/* functions in hor.c */
void mps_horner(mps_status* s, mpc_t y, int *dprec, int *iprec, int deg, int i);
void mps_refine(mps_status* s, int i, long prec);
void mps_ruroutroot(mps_status* s, mpc_t root, char status, long prec, long out_prec);


/*
 * End of extern "C" {
 *   ...
 * }
 */
#ifdef __cplusplus
}
#endif

#endif /* ndef MPS_RURSOLVE_H */