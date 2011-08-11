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

#include <mps/core.h>
#include <mps/secular.h>

static long int data_prec_max = 0;

/***********************************************************
 *           SUBROUTINE MP_SET_PREC                        *
 ***********************************************************
  Globally set the current precision of mp variables
 **********************************************************/
void
mps_mp_set_prec (mps_status * s, long int prec)
{
  s->mpwp = prec;
  mpf_set_default_prec (prec);
  rdpe_set_2dl (s->mp_epsilon, 1.0, -prec + 1);
}

/********************************************************
 *      SUBROUTINE ALLOCATE_DATA                        *
 *******************************************************/
void
mps_allocate_data (mps_status * s)
{
  int i;

  if (s->DOLOG)
    fprintf (s->logstr, "Allocating data\n");

  s->clust = int_valloc (s->deg);
  s->punt = int_valloc (s->deg + 1);
  s->clust_detached = int_valloc (s->deg);

  s->again = mps_boolean_valloc (s->deg);

  s->status = (char (*)[3]) char_valloc (3 * s->deg);

  s->order = int_valloc (s->deg);
  s->rootwp = long_valloc (s->deg);

  s->fap = double_valloc (s->deg + 1);
  s->dap = rdpe_valloc (s->deg + 1);

  s->frad = double_valloc (s->deg);
  s->froot = cplx_valloc (s->deg);
  s->drad = rdpe_valloc (s->deg);
  s->droot = cdpe_valloc (s->deg);

  s->mroot = mpc_valloc (s->deg);
  for (i = 0; i < s->deg; i++)
    mpc_init2 (s->mroot[i], 0);

  s->fppc = cplx_valloc (s->deg + 1);
  s->fppc1 = cplx_valloc (s->deg + 1);

  s->mfpc1 = mpc_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    mpc_init2 (s->mfpc1[i], 0);

  s->mfppc = mpc_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    mpc_init2 (s->mfppc[i], 0);

  s->mfppc1 = mpc_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    mpc_init2 (s->mfppc1[i], 0);

  s->mfpc2 = mpc_valloc ((s->deg + 1) * s->n_threads);
  for (i = 0; i < (s->deg + 1) * s->n_threads; i++)
    mpc_init2 (s->mfpc2[i], 0);

  /* temporary vectors */
  s->spar1 = mps_boolean_valloc (s->deg + 2);
  s->spar2 = mps_boolean_valloc ((s->deg + 2) * s->n_threads);
  s->h = mps_boolean_valloc (s->deg + 2);
  s->again_old = mps_boolean_valloc (s->deg);

  s->oldpunt = int_valloc (s->deg + 1);
  s->clust_aux = int_valloc (s->deg + 1);
  s->punt_aux = int_valloc (s->deg + 1);
  s->punt_out = int_valloc (s->deg + 1);
  s->clust_out = int_valloc (s->deg + 1);

  s->fap1 = double_valloc (s->deg + 1);
  s->fap2 = double_valloc (s->deg + 1);

  s->dap1 = rdpe_valloc (s->deg + 1);
  s->dap2 = rdpe_valloc ((s->deg + 1) * s->n_threads);
  s->dpc1 = cdpe_valloc (s->deg + 1);
  s->dpc2 = cdpe_valloc (s->deg + 1);

  s->fradii = double_valloc (s->deg + 1);
  s->partitioning = int_valloc (s->deg + 2);
  s->dradii = rdpe_valloc (s->deg + 1);
}

/***********************************************************
 *           SUBROUTINE RAISE_DATA                         *
 ***********************************************************
 raise precision performing a real computation of the data
 **********************************************************/
long int
mps_raise_data (mps_status * s, long int prec)
{
  int i, k;

  /* raise the precision of  mroot */
  for (k = 0; k < s->n; k++)
    mpc_set_prec (s->mroot[k], prec);

  if (s->data_type[0] != 'u')
    {
      /* raise the precision of  mfpc */
      for (k = 0; k < s->n + 1; k++)
	if (s->data_type[0] != 's' || s->spar[k])
	  mpc_set_prec (s->mfpc[k], prec);

      for (i = 0; i <= s->n; i++)
	if (s->data_type[0] != 's' || s->spar[i])
	  {
	    switch (s->data_type[1])
	      {
	      case 'r':	/* real */
		/* the real case should be adjusted later on */
		switch (s->data_type[2])
		  {
		  case 'i':	/* integer */
		    mpf_set_z (mpc_Re (s->mfpc[i]), s->mip_r[i]);
		    mpf_set_ui (mpc_Im (s->mfpc[i]), 0);
		    break;
		  case 'q':	/* rational */
		    mpf_set_q (mpc_Re (s->mfpc[i]), s->mqp_r[i]);
		    mpf_set_ui (mpc_Im (s->mfpc[i]), 0);
		    /* GMP 2.0.2 bug begin */
		    if (mpf_sgn (mpc_Re (s->mfpc[i])) !=
			mpq_sgn (s->mqp_r[i]))
		      mpf_neg (mpc_Re (s->mfpc[i]), mpc_Re (s->mfpc[i]));
		    /* GMP bug end */
		    break;
		  case 'b':	/* big float */
		    break;	/* nothing to do */
		  case 'f':	/* float */
		    /* real case
		       mpc_set_d(mfpc[i], fpr[i], 0.0);
		     */
		    break;
		  default:
		    mps_error (s, 1, "Mistake in goal");
		    break;
		  }
		break;
	      case 'c':	/* complex */
		switch (s->data_type[2])
		  {
		  case 'i':	/* integer */
		    mpc_set_z (s->mfpc[i], s->mip_r[i], s->mip_i[i]);
		    break;
		  case 'q':	/* rational */
		    mpc_set_q (s->mfpc[i], s->mqp_r[i], s->mqp_i[i]);
		    /* GMP 2.0.2 bug begin */
		    if (mpf_sgn (mpc_Re (s->mfpc[i])) !=
			mpq_sgn (s->mqp_r[i]))
		      mpf_neg (mpc_Re (s->mfpc[i]), mpc_Re (s->mfpc[i]));
		    if (mpf_sgn (mpc_Im (s->mfpc[i])) !=
			mpq_sgn (s->mqp_i[i]))
		      mpf_neg (mpc_Im (s->mfpc[i]), mpc_Im (s->mfpc[i]));
		    /* GMP bug end */
		    break;
		  case 'b':	/* big float */
		    break;	/* nothing to do */
		  case 'f':	/* float */
		    /* real case
		       mpc_set_cplx(mfpc[i], fpc[i]);
		     */
		    break;
		  default:
		    mps_error (s, 1, "Mistake in goal");
		    break;
		  }
	      }
	  }
    }

  /* Raise the precision of p' */
  if (s->data_type[0] == 's')
    for (k = 0; k < s->n; k++)
      if (s->spar[k + 1])
	{
	  mpc_set_prec (s->mfppc[k], prec);
	  mpc_mul_ui (s->mfppc[k], s->mfpc[k + 1], k + 1);
	}

  /* raise the precision of auxiliary variables */
  for (k = 0; k < s->n + 1; k++)
    {
      mpc_set_prec (s->mfpc1[k], prec);
      mpc_set_prec (s->mfppc1[k], prec);
    }

  if (s->data_type[0] == 's')
    for (k = 0; k < (s->n + 1) * s->n_threads; k++)
      mpc_set_prec (s->mfpc2[k], prec);

  return mpc_get_prec (s->mroot[0]);
}

/***********************************************************
 *           SUBROUTINE RAISE_DATA_RAW                     *
 ***********************************************************
 modify the raw precision of mp variables
 ***********************************************************/
void
mps_raise_data_raw (mps_status * s, long int prec)
{
  int k;

  /* raise the precision of  mroot */
  for (k = 0; k < s->n; k++)
    mpc_set_prec_raw (s->mroot[k], prec);

  /* raise the precision of  mfpc */
  if (s->data_type[0] != 'u')
    for (k = 0; k < s->n + 1; k++)
      mpc_set_prec_raw (s->mfpc[k], prec);

  /* Raise the precision of sparse vectors */
  if (s->data_type[0] == 's')
    for (k = 0; k < s->n; k++)
      if (s->spar[k + 1])
	mpc_set_prec_raw (s->mfppc[k], prec);

  /* raise the precision of auxiliary variables */
  for (k = 0; k < s->n + 1; k++)
    {
      mpc_set_prec_raw (s->mfpc1[k], prec);
      mpc_set_prec_raw (s->mfppc1[k], prec);
    }

  if (s->data_type[0] == 's')
    for (k = 0; k < (s->n + 1) * s->n_threads; k++)
      mpc_set_prec_raw (s->mfpc2[k], prec);
}

/***********************************************************
 *           SUBROUTINE PREPARE_DATA                       *
 ***********************************************************
 Compute the mp_complex values of the coefficients of p(x)
 with the  current precision of mpwds words, given the
 rational or integer coefficients.
 ***********************************************************/
void
mps_prepare_data (mps_status * s, long int prec)
{
  MPS_DEBUG_THIS_CALL

  MPS_DEBUG(s, "Increasing working precision to %ld bits", prec);

  if (prec > data_prec_max)
    {
      if (data_prec_max)
	mps_raise_data_raw (s, data_prec_max);
      data_prec_max = mps_raise_data (s, prec);
    }
  else
  {
      /* Check if the algorithm is Standard MPSolve or the secular
       * equation version */
      if (s->mpsolve_ptr == MPS_MPSOLVE_PTR (mps_standard_mpsolve))
          mps_raise_data_raw (s, prec);
      else
          mps_secular_raise_precision (s, prec);
  }
}

/***********************************************************
 *           SUBROUTINE RESTORE_DATA                       *
 ***********************************************************
 Resets the data to the highest used precision
 ***********************************************************/
void
mps_restore_data (mps_status * s)
{
  if (s->DOLOG)
    fprintf (s->logstr, "Restore data to %ld bits\n", data_prec_max);

  if (data_prec_max)
    mps_raise_data_raw (s, data_prec_max);
}

/********************************************************
 *      SUBROUTINE FREE_DATA                            *
 *******************************************************/
void
mps_free_data (mps_status * s)
{
  int i;

  if (s->DOLOG)
    fprintf (s->logstr, "Unallocating data...\n");

  free (s->clust);
  free (s->punt);
  free (s->clust_detached);
  free (s->again);
  free (s->status);
  free (s->rootwp);
  free (s->order);

  free (s->fap);
  rdpe_vfree (s->dap);

  free (s->frad);
  rdpe_vfree (s->drad);

  cplx_vfree (s->froot);
  cdpe_vfree (s->droot);
  for (i = 0; i < s->deg; i++)
    mpc_clear (s->mroot[i]);
  free (s->mroot);

  for (i = 0; i <= s->deg; i++)
    mpc_clear (s->mfpc1[i]);
  mpc_vfree (s->mfpc1);

  cplx_vfree (s->fppc);
  cplx_vfree (s->fppc1);
  for (i = 0; i <= s->deg; i++)
    {
      mpc_clear (s->mfppc[i]);
      mpc_clear (s->mfppc1[i]);
    }

  for (i = 0; i < (s->deg + 1) * s->n_threads; i++)
    {
      mpc_clear (s->mfpc2[i]);
    }

  free (s->mfppc);
  free (s->mfppc1);
  free (s->mfpc2);


  /* free temporary vectors */
  free (s->spar1);
  mps_boolean_vfree (s->spar2);
  free (s->h);
  free (s->again_old);

  free (s->oldpunt);
  free (s->clust_aux);
  free (s->punt_aux);
  free (s->punt_out);
  free (s->clust_out);

  free (s->fap1);
  free (s->fap2);

  rdpe_vfree (s->dap1);
  rdpe_vfree (s->dap2);
  cdpe_vfree (s->dpc1);
  cdpe_vfree (s->dpc2);

  free (s->partitioning);
  free (s->fradii);
  rdpe_vfree (s->dradii);

  if (s->DOLOG)
    fprintf (s->logstr, "...temporaries...\n");

  mptemp_clear ();

  if (s->DOLOG)
    fprintf (s->logstr, "...done\n");
}
