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

#include "mps.h"

static long int data_prec_max = 0;

/***********************************************************
 *           SUBROUTINE MP_SET_PREC                        *
 ***********************************************************
  Globally set the current precision of mp variables
 **********************************************************/
void
mp_set_prec(long int prec)
{
  mpwp = prec;
  mpf_set_default_prec(prec);
  rdpe_set_2dl(mp_epsilon, 1.0, -prec + 1);
}

/********************************************************
 *      SUBROUTINE ALLOCATE_DATA                        *
 *******************************************************/
void
allocate_data(void)
{
  int i;

  if (DOLOG)
    fprintf(logstr, "Allocating data\n");

  clust = int_valloc(deg);
  punt = int_valloc(deg + 1);
  again = boolean_valloc(deg);
  status = (char (*)[3]) char_valloc(3 * deg);
  
  order = int_valloc(deg);
  rootwp = long_valloc(deg);

  fap = double_valloc(deg + 1);
  dap = rdpe_valloc(deg + 1);

  frad = double_valloc(deg);
  froot = cplx_valloc(deg);
  drad = rdpe_valloc(deg);
  droot = cdpe_valloc(deg);

  mroot = mpc_valloc(deg);
  for (i = 0; i < deg; i++)
    mpc_init2(mroot[i], 0);

  fppc = cplx_valloc(deg + 1);
  fppc1 = cplx_valloc(deg + 1);

  mfpc1 = mpc_valloc(deg + 1);
  for (i = 0; i <= deg; i++)
    mpc_init2(mfpc1[i], 0);

  mfppc = mpc_valloc(deg + 1);
  mfppc1 = mpc_valloc(deg + 1);
  mfpc2 = mpc_valloc(deg + 1);
  for (i = 0; i <= deg; i++) {
    mpc_init2(mfppc[i], 0);
    mpc_init2(mfppc1[i], 0);
    mpc_init2(mfpc2[i], 0);
  }

  /* temporary vectors */
  spar1 = boolean_valloc(deg + 2);
  spar2 = boolean_valloc(deg + 2);
  h = boolean_valloc(deg + 2);
  again_old = boolean_valloc(deg);

  oldpunt = int_valloc(deg + 1);
  clust_aux = int_valloc(deg + 1);
  punt_aux = int_valloc(deg + 1);
  punt_out = int_valloc(deg + 1);
  clust_out = int_valloc(deg + 1);

  fap1 = double_valloc(deg + 1);
  fap2 = double_valloc(deg + 1);

  dap1 = rdpe_valloc(deg + 1);
  dap2 = rdpe_valloc(deg + 1);
  dpc1 = cdpe_valloc(deg + 1);
  dpc2 = cdpe_valloc(deg + 1);
}

/***********************************************************
 *           SUBROUTINE RAISE_DATA                         *
 ***********************************************************
 raise precision performing a real computation of the data
 **********************************************************/
long int
raise_data(long int prec)
{
  int i, k;

  /* raise the precision of  mroot */
  for (k = 0; k < n; k++)
    mpc_set_prec(mroot[k], prec);

  if (data_type[0] != 'u') {
    /* raise the precision of  mfpc */
    for (k = 0; k < n + 1; k++)
      if (data_type[0] != 's' || spar[k])
	mpc_set_prec(mfpc[k], prec);

    for (i = 0; i <= n; i++)
      if (data_type[0] != 's' || spar[i]) {
	switch (data_type[1]) {
	case 'r':		/* real */
	  /* the real case should be adjusted later on */
	  switch (data_type[2]) {
	  case 'i':		/* integer */
	    mpf_set_z(mpc_Re(mfpc[i]), mip_r[i]);
	    mpf_set_ui(mpc_Im(mfpc[i]), 0);
	    break;
	  case 'q':		/* rational */
	    mpf_set_q(mpc_Re(mfpc[i]), mqp_r[i]);
	    mpf_set_ui(mpc_Im(mfpc[i]), 0);
	    /* GMP 2.0.2 bug begin */
	    if (mpf_sgn(mpc_Re(mfpc[i])) != mpq_sgn(mqp_r[i]))
	      mpf_neg(mpc_Re(mfpc[i]), mpc_Re(mfpc[i]));
	    /* GMP bug end */
	    break;
	  case 'b':		/* big float */
	    break;		/* nothing to do */
	  case 'f':		/* float */
	    /* real case
	    mpc_set_d(mfpc[i], fpr[i], 0.0);
	    */
	    break;
	  default:
	    error(1, "Mistake in goal");
	    break;
	  }
	  break;
	case 'c':		/* complex */
	  switch (data_type[2]) {
	  case 'i':		/* integer */
	    mpc_set_z(mfpc[i], mip_r[i], mip_i[i]);
	    break;
	  case 'q':		/* rational */
	    mpc_set_q(mfpc[i], mqp_r[i], mqp_i[i]);
	    /* GMP 2.0.2 bug begin */
	    if (mpf_sgn(mpc_Re(mfpc[i])) != mpq_sgn(mqp_r[i]))
	      mpf_neg(mpc_Re(mfpc[i]), mpc_Re(mfpc[i]));
	    if (mpf_sgn(mpc_Im(mfpc[i])) != mpq_sgn(mqp_i[i]))
	      mpf_neg(mpc_Im(mfpc[i]), mpc_Im(mfpc[i]));
	    /* GMP bug end */
	    break;
	  case 'b':		/* big float */
	    break;		/* nothing to do */
	  case 'f':		/* float */
	    /* real case
	    mpc_set_cplx(mfpc[i], fpc[i]);
	    */
	    break;
	  default:
	    error(1, "Mistake in goal");
	    break;
	  }
	}
      }
  }
  
  /* Raise the precision of p' */
  if (data_type[0] == 's')
    for (k = 0; k < n; k++)
      if (spar[k + 1]) {
	mpc_set_prec(mfppc[k], prec);
	mpc_mul_ui(mfppc[k], mfpc[k + 1], k + 1);
      }

  /* raise the precision of auxiliary variables */
  for (k = 0; k < n + 1; k++) {
    mpc_set_prec(mfpc1[k], prec);
    mpc_set_prec(mfppc1[k], prec);
  }
  
  if (data_type[0] == 's')
    for (k = 0; k < n + 1; k++)
      mpc_set_prec(mfpc2[k], prec);

  return mpc_get_prec(mroot[0]);
}

/***********************************************************
 *           SUBROUTINE RAISE_DATA_RAW                     *
 ***********************************************************
 modify the raw precision of mp variables
 ***********************************************************/
void
raise_data_raw(long int prec)
{
  int k;

  /* raise the precision of  mroot */
  for (k = 0; k < n; k++)
    mpc_set_prec_raw(mroot[k], prec);

  /* raise the precision of  mfpc */
  if (data_type[0] != 'u')
    for (k = 0; k < n + 1; k++)
      mpc_set_prec_raw(mfpc[k], prec);

  /* Raise the precision of sparse vectors */
  if (data_type[0] == 's')
    for (k = 0; k < n; k++)
      if (spar[k + 1])
	mpc_set_prec_raw(mfppc[k], prec);

  /* raise the precision of auxiliary variables */
  for (k = 0; k < n + 1; k++) {
    mpc_set_prec_raw(mfpc1[k], prec);
    mpc_set_prec_raw(mfppc1[k], prec);
  }

  if (data_type[0] == 's')
    for (k = 0; k < n + 1; k++)
      mpc_set_prec_raw(mfpc2[k], prec);
}

/***********************************************************
 *           SUBROUTINE PREPARE_DATA                       *
 ***********************************************************
 Compute the mp_complex values of the coefficients of p(x)
 with the  current precision of mpwds words, given the
 rational or integer coefficients.
 ***********************************************************/
void
prepare_data(long int prec)
{
  if (DOLOG)
    fprintf(logstr, "Prepare data:  working precision =%ld bits\n", prec);

  if (prec > data_prec_max) {
    if (data_prec_max)
      raise_data_raw(data_prec_max);
    data_prec_max = raise_data(prec);
  } else
    raise_data_raw(prec);
}

/***********************************************************
 *           SUBROUTINE RESTORE_DATA                       *
 ***********************************************************
 Resets the data to the highest used precision
 ***********************************************************/
void
restore_data(void)
{
  if (DOLOG)
    fprintf(logstr, "Restore data to %ld bits\n", data_prec_max);

  if (data_prec_max) 
    raise_data_raw(data_prec_max);
}

/********************************************************
 *      SUBROUTINE FREE_DATA                            *
 *******************************************************/
void
free_data(void)
{
  int i;

  if (DOLOG)
    fprintf(logstr, "Unallocating data...\n");

  free(clust);
  free(punt);
  free(again);
  free(status);
  free(rootwp);
  free(order);

  free(fap);
  rdpe_vfree(dap);

  free(frad);
  rdpe_vfree(drad);

  cplx_vfree(froot);
  cdpe_vfree(droot);
  for (i = 0; i < deg; i++)
    mpc_clear(mroot[i]);
  free(mroot);

  for (i = 0; i <= deg; i++)
    mpc_clear(mfpc1[i]);
  free(mfpc1);

  cplx_vfree(fppc1);
  for (i = 0; i <= deg; i++) {
    mpc_clear(mfppc[i]);
    mpc_clear(mfppc1[i]);
    mpc_clear(mfpc2[i]);
  }
  free(mfppc);
  free(mfppc1);
  free(mfpc2);


  /* free temporary vectors */
  free(spar1);
  free(spar2);
  free(h);
  free(again_old);

  free(oldpunt);
  free(clust_aux);
  free(punt_aux);
  free(punt_out);
  free(clust_out);

  free(fap1);
  free(fap2);

  rdpe_vfree(dap1);
  rdpe_vfree(dap2);
  cdpe_vfree(dpc1);
  cdpe_vfree(dpc2);

  if (DOLOG)
    fprintf(logstr, "...temporaries...\n");

  mptemp_clear();

  if (DOLOG)
    fprintf(logstr, "...done\n");
}
