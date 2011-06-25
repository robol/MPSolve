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

#ifndef MPS_POLY_H
#define MPS_POLY_H

#ifdef __cplusplus
extern "C" {
#endif


#include <mps/core.h>

/* functions for a simplified I/O with the mpsolve library */


/**
 * @brief struct for the polynomial data type.
 */
typedef struct {
  int deg;		/* starting polynomial degree */
  char data_type[3];	/* polynomial data type */
  long int prec_in;	/* number of digits of input precision */
  int n;		/* degree */
  mps_boolean *spar;	/* sparsity structure of the polynomial */
  double *fpr;		/* standard real coefficients */
  cplx_t *fpc;		/* standard complex coefficients */
  rdpe_t *dpr;		/* dpe real coefficients */
  cdpe_t *dpc;		/* dpe complex coefficients */
  mpz_t *mip_r;		/* real part of integer input coefs */
  mpz_t *mip_i;		/* imaginary part of integer input coefs */
  mpq_t *mqp_r;		/* real part of rational input coeff. */
  mpq_t *mqp_i;		/* imaginary part of rational input coefs */
  mpf_t *mfpr;		/* multiprecision real coefficients */
  mpc_t *mfpc;		/* multiprecision complex coefficients */
} __mpspoly_struct;


/**
 * @brief pointer to polynomial data type __mpspoly_struct.
 */
typedef __mpspoly_struct mpspoly_t[1];

/* solution data type */
typedef struct {
  mps_phase lastphase;	/* store last computed phase */
  int *count;		/* count roots: [inside, outside, uncertain] */
  int zero_roots;	/* number of roots = 0 */
  char (*status)[3];	/* status of each approximation */
  cplx_t *froot;	/* root approx. as standard complex numbers */
  cdpe_t *droot;	/* root approximations as complex dpe numbers */
  mpc_t *mroot;		/* root approximations as complex mp numbers */
  double *frad;		/* radii of the incl. disks as real numbers */
  rdpe_t *drad;		/* radii of the incl. disks as rdpe_t numbers */
} __mpsroots_struct;

typedef __mpsroots_struct mpsroots_t[1];

/* functions */
void mps_allocate_poly(mps_status* s, mpspoly_t p);
void mps_read_poly(mps_status* s, FILE *instr, mpspoly_t p);
void mps_validate_poly(mps_status* s, mpspoly_t p, int num_coeff);
void mps_set_poly(mps_status* s, mpspoly_t p);
void mps_update_poly(mps_status* s, mpspoly_t p);
void mps_free_poly(mps_status* s, mpspoly_t p);
void mps_get_roots(mps_status* s, mpsroots_t r);

/*
 * End of extern "C" {
 *   ...
 * }
 */
#ifdef __cplusplus
}
#endif

#endif /* ndef MPS_POLY_H */
