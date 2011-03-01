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

#include "mps.h"

extern mpz_t *mpdemo;	/* imaginary part of the integer input coeff. */

/* functions in main.c */
void rursolve(void);

/* functions in hor.c */
void horner(mpc_t y, int *dprec, int *iprec, int deg, int i);
void refine(int i, long prec);
void ruroutroot(mpc_t root, char status, long prec, long out_prec);
