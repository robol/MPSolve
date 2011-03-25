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

/***********************************************************
 *           SUBROUTINE MPSolve                            *
 ***********************************************************
 The program is divided into many parts
 - Check the correctness of data, scale coefficients if
   needed, and select cases: the variable which_case is
   'f' or 'd'  according to float or dpe case.
 - Call msolve or dsolve according to the value of which_case.
 - Allocate MP variables mfpc, mroot, drad (if needed).
 - Start MPsolve loop
     - prepare data according to the current precision
       and to the data_type (density/sparsity/user)  
     - Call msolve with the current precision
 - check for termination
************************************************************/
void
mpsolve(void)
{
  int i, nzc;
  char which_case;
  boolean d_after_f, computed, over_max;

  /* == 1 ==  Setup variables, i.e. copy coefficients
   into dpr, dpc and similar. */
  setup();

  lastphase = no_phase;
  computed = false;
  over_max = false;
  
  /* == 2 ==  Resume from pre-computed roots */
  if (resume) /* to complete */
    error(1, "Resume not supported yet");
  
  /* == 3 ==  Check data and get starting phase */
  if (skip_float)
    which_case = 'd';
  else
    which_case = 'f';
  d_after_f = false;
  
  check_data(&which_case);
  if (DOLOG)
    fprintf(logstr, "Which_case = %c, skip_float= %d\n",
	    which_case, skip_float);

  /* == 4 ==  Float phase */
  if (which_case == 'f') {
    if (DOLOG)
      fprintf(logstr, "Float phase ...\n");
    fsolve(&d_after_f);
    lastphase = float_phase;

    if (DOLOG)
      dump(logstr);

    computed = check_stop();
    if (computed && goal[0] != 'a')
      goto exit_sub;  /* stop for COUNT and ISOLATE goals */
  }

  /* == 5 ==  DPE phase */
  if (which_case == 'd' || d_after_f) {	/* DPE phase */
    if (DOLOG)
      fprintf(logstr, "DPE phase ...\n");
    if (d_after_f)
      for (i = 0; i < n; i++) {
	rdpe_set_d(drad[i], frad[i]);
	cdpe_set_x(droot[i], froot[i]);
      }
    dsolve(d_after_f);
    lastphase = dpe_phase;

    if (DOLOG)
      dump(logstr);

    computed = check_stop();
    if (computed && goal[0] != 'a')
      goto exit_sub;
  }

  /* == 6 ==   Allocate MP variables mfpc, mroot, drad, mfppc, mfppc1
   * (the real input case is not implemented yet ) */
  if (DOLOG)
    fprintf(logstr, "MP phase ...\n");

  /* ==== 6.1 initialize mp variables */
  mp_set_prec(2 * DBL_MANT_DIG);

  /* Prepare data according to the current working precision */
  prepare_data(mpwp);

  /* ==== 6.2 set initial values for mp variables */
  for (i = 0; i < n; i++)
    if (which_case == 'd' || d_after_f)
      mpc_set_cdpe(mroot[i], droot[i]);
    else {
      mpc_set_cplx(mroot[i], froot[i]);
      rdpe_set_d(drad[i], frad[i]);
    }
  if (computed && goal[0] == 'a')
    goto exit_sub;
  
  /* == 7 ==  Start MPsolve loop */
  mpwp = DBL_MANT_DIG;
  while (!computed && mpwp < mpwp_max && (prec_in == 0 || mpwp < prec_in)) {

    mpwp *= 2;

    if (prec_in != 0 && mpwp > prec_in)
      mpwp = prec_in + (int) (log(4.0 * n) / LOG2);

    if (mpwp > mpwp_max) {
      mpwp = mpwp_max;
      over_max = true;
    }

    if (DOLOG)
      fprintf(logstr, "MAIN: mp_loop: mpwp=%ld\n", mpwp);

    /* == 7.1 ==   prepare data according to the current precision */
    mp_set_prec(mpwp);
    prepare_data(mpwp);

    /* == 7.2 ==   Call msolve with the current precision */
    if (DOLOG)
      fprintf(logstr, "MAIN: now call msolve nclust=%d\n", nclust);
    msolve();
    lastphase = mp_phase;

    /* if (DOLOG) dump(logstr); */
    
    if (DOLOG) {  /* count isolated zeros */
      nzc = 0;
      for (i = 0; i < n; i++) {
	if (status[i][0] == 'i' || status[i][0] == 'a')
	  nzc++;
      }     
      fprintf(logstr, "MAIN: isolated %d roots\n", nzc);
      fprintf(logstr, "MAIN: after msolve check stop\n");
    }
    
    /* == 7.3 ==  Check the stop condition */
    computed = check_stop();
    mmodify();

    /* == 7.4 ==  reset the status vector */
    for (i = 0; i < n; i++)
      if (status[i][0] == 'C')
	status[i][0] = 'c';
  }

  /* == 8 ==  Check for termination */
  if (!computed) {
    if (over_max)
      error(1, "Reached the maximum working precision");
    else
      warn("Reached the input precision");
  }
  
 exit_sub:

  /* == 9 ==  Check inclusion disks */
  if (computed && nclust<n)
    if (!inclusion())
      error(1, "Unable to compute inclusion disks");
  
  /* == 10 ==  Refine roots */
  if (computed && !over_max && goal[0] == 'a') {
    lastphase = mp_phase;
    improve();
  }
  
  /* == 11 ==  Restore to highest used precision */
  if (lastphase == mp_phase)
    restore_data();
}

/***********************************************************
 *           SUBROUTINE SETUP                              *
 ***********************************************************
 Setup vectors and variables
 ***********************************************************/
void
setup(void)
{
  int i;
  tmpf_t mptemp;
  tmpc_t mptempc;
  
  if (DOLOG) {
    fprintf(logstr, "Goal      = %5s\n", goal);
    fprintf(logstr, "Data type = %3s\n", data_type);
    fprintf(logstr, "Degree    = %d\n", n);
    fprintf(logstr, "Input prec.  = %ld digits\n", (long) (prec_in * LOG10_2));
    fprintf(logstr, "Output prec. = %ld digits\n", (long) (prec_out * LOG10_2));
  }

  /* setup temporary vectors */
  if (data_type[0] == 's')
    for (i = 0; i <= n; i++) {
      fap[i] = 0.0;
      fpr[i] = 0.0;
      rdpe_set(dap[i], rdpe_zero);
      cplx_set(fpc[i], cplx_zero);
      rdpe_set(dpr[i], rdpe_zero);
      cdpe_set(dpc[i], cdpe_zero);      
     }
  
  /* setup status and clusters so that there is only one cluster
  *  containing all the roots */
  for (i = 0; i < n; i++) {
    /* Set the i-th root as a clustered root */
    status[i][0] = 'c';
    /* Set the i-th root as 'uncertain', regarding the
     * inclusion in R or in iR. */
    status[i][1] = 'w';
    /* Set the i-th root as 'uncertain' (regarding the inclusion
     * in the search set). */
    status[i][2] = 'u';
    clust[i] = i;
  }

  /* Indexes of the first (and only) cluster start from
   * 0 and reach n */
  punt[0] = 0;
  punt[1] = n;
  
  /* set input and output epsilon */
  rdpe_set_2dl(eps_in, 1.0, 1 - prec_in);
  rdpe_set_2dl(eps_out, 1.0, 1 - prec_out);

  /* precision of each root */
  for (i = 0; i < n; i++)
    rootwp[i] = 53;
   
  /* output order info */
  for (i = 0; i < n; i++)
    order[i] = i;

  /* reset root counts */
  count[0] = count[1] = count[2] = 0;

  /* compute DPE approximations */
  if (data_type[0] == 'u')
    return; /* nothing to do */
    
  /* init temporary mp variables */
  tmpf_init2(mptemp, DBL_MANT_DIG);
  tmpc_init2(mptempc, DBL_MANT_DIG);

  /* main loop */
  skip_float = false;
  for (i = 0; i <= n; i++) {

    if (data_type[0] == 's' && !spar[i])
      continue;

    switch (data_type[1]) {	/* switch 1 */

    case 'r':			/* Real */

      switch (data_type[2]) {	/* switch 2 */

      case 'i':		/* Real - Integer Coefs */
	mpf_set_z(mptemp, mip_r[i]);
	mpf_get_rdpe(dpr[i], mptemp);
	break;

      case 'q':		/* Real - Rational Coefs */
	mpf_set_q(mptemp, mqp_r[i]);
	mpf_get_rdpe(dpr[i], mptemp);
	/*#G GMP 2.0.2 bug begin */
	if (rdpe_sgn(dpr[i]) != mpq_sgn(mqp_r[i]))
	  rdpe_neg_eq(dpr[i]);
	/*#G GMP bug end */
	break;

      case 'f':		/* Real - Big/Float Coefs */
	mpf_get_rdpe(dpr[i], mfpr[i]);
	break;

      }				/* switch 2 */

      cdpe_set_e(dpc[i], dpr[i], rdpe_zero);
      
      /* compute dap[i] and check for float phase */
      rdpe_abs(dap[i], dpr[i]);
      rdpe_abs(dap[i], dpr[i]);
      if (rdpe_gt(dap[i], rdpe_maxd) || rdpe_lt(dap[i], rdpe_mind))
	skip_float = true;
	
      break;			/* switch 1 */

    case 'c':			/* Complex */

      switch (data_type[2]) {	/* switch 3 */

      case 'i':		/* Complex - Integer Coefs */
	mpc_set_z(mptempc, mip_r[i], mip_i[i]);
	mpc_get_cdpe(dpc[i], mptempc);
	break;

      case 'q':		/* Complex - Rational Coefs */
	mpc_set_q(mptempc, mqp_r[i], mqp_i[i]);
	mpc_get_cdpe(dpc[i], mptempc);
	/*#G GMP 2.0.2 bug begin */
	if (rdpe_sgn(cdpe_Re(dpc[i])) != mpq_sgn(mqp_r[i]))
	  rdpe_neg_eq(cdpe_Re(dpc[i]));
	if (rdpe_sgn(cdpe_Im(dpc[i])) != mpq_sgn(mqp_i[i]))
	  rdpe_neg_eq(cdpe_Im(dpc[i]));
	/*#G GMP bug end */
	break;

      case 'f':		/* Complex - Big/Float Coefs */
	mpc_get_cdpe(dpc[i], mfpc[i]);
	break;

      }				/* switch 3 */

      /* compute dap[i] */
      cdpe_mod(dap[i], dpc[i]);
      if (rdpe_gt(dap[i], rdpe_maxd) || rdpe_lt(dap[i], rdpe_mind))
	skip_float = true;

      break;

    }				/* switch 1 */
  }				/* for */

  /* free temporary mp variables */
  tmpf_clear(mptemp);
  tmpc_clear(mptempc);

  /* adjust input data type */
  if (data_type[2] == 'f' && skip_float) 
    data_type[2] = 'b';

  /* prepare floating point vectors */
  if (!skip_float) 
    for (i = 0; i <= n; i++) {
      if (data_type[0] == 's' || !spar[i])
	continue;
      if (data_type[1] == 'r') {
	fpr[i] = rdpe_get_d(dpr[i]);
	fap[i] = fabs(fpr[i]);
	cplx_set_d(fpc[i], fpr[i], 0.0);
      } else {
	cdpe_get_x(fpc[i], dpc[i]);
	fap[i] = cplx_mod(fpc[i]);
      }
    }
}

/*********************************************************
*                   PROCEDURE CHECK_DATA
**********************************************************/
void
check_data(char *which_case)
{
  rdpe_t min_coeff, max_coeff, tmp;
  cdpe_t ctmp;
  int i;

  /* case of user-defined polynomial */
  if (data_type[0] == 'u') {
    if (goal[2] == 'm')
      error(1, "Multiplicity detection not yet implemented for user polynomial");
    if (goal[3] != 'n')
      error(1, "Real/imaginary detection not yet implemented for user polynomial");
    *which_case = 'd';
    return;
  }

  /* Check consistency of input */
  if (rdpe_eq(dap[n], rdpe_zero)) {
    warn("The leading coefficient is zero");
    do
      n--;
    while (rdpe_eq(dap[n], rdpe_zero));
  }

  /* count number of zero roots (undeflated input polynomial) */
  zero_roots = 0;
  while (rdpe_eq(dap[zero_roots], rdpe_zero))
    zero_roots++;
  /* shift down input vectors */
  if (zero_roots) {
    for (i = 0; i <= n - zero_roots; i++) {
      rdpe_set(dap[i], dap[i + zero_roots]);
      fap[i] = fap[i + zero_roots];
      fpr[i] = fpr[i + zero_roots];
      cplx_set(fpc[i], fpc[i + zero_roots]);
      rdpe_set(dpr[i], dpr[i + zero_roots]);
      cdpe_set(dpc[i], dpc[i + zero_roots]);
      mpf_set(mfpr[i], mfpr[i + zero_roots]);
      mpc_set(mfpc[i], mfpc[i + zero_roots]);
      if (i < n - zero_roots)
	mpc_set(mfppc[i], mfppc[i + zero_roots]);
      mpz_set(mip_r[i], mip_r[i + zero_roots]);
      mpz_set(mip_i[i], mip_i[i + zero_roots]);
      mpq_set(mqp_r[i], mqp_r[i + zero_roots]);
      mpq_set(mqp_i[i], mqp_i[i + zero_roots]);
      spar[i] = spar[i + zero_roots];
    }
    n = n - zero_roots;
  }

  /* Compute min_coeff */
  if (rdpe_lt(dap[0], dap[n]))
    rdpe_set(min_coeff, dap[0]);
  else
    rdpe_set(min_coeff, dap[n]);

  /* Compute max_coeff and its logarithm */
  rdpe_set(max_coeff, dap[0]);
  for (i = 1; i <= n; i++)
    if (rdpe_lt(max_coeff, dap[i]))
      rdpe_set(max_coeff, dap[i]);
  lmax_coeff = rdpe_log(max_coeff);

  /*  Multiplicity and sep */
  if (goal[2] == 'm')
    switch (data_type[2]) {
    case 'i':
      compute_sep();
      break;
    case 'q':
      warn("The multiplicity option has not been yet implemented");
      sep = 0.0;
      break;
    default:
      warn("The input polynomial has neither integer nor rational");
      warn(" coefficients: unable to compute multiplicities");
      sep = 0.0;
      break;
    }

  /* Real/Imaginary detection */
  if (goal[3] != 'n' || goal[1] == 'R' || goal[1] == 'I')
    switch (data_type[2]) {
    case 'i':
      compute_sep();
      break;
    case 'q':
      warn("The  real/imaginary option has not been yet implemented");
      sep = 0.0;
      break;
    default:
      warn("The input polynomial has neither integer nor rational");
      warn(" coefficients: unable to perform real/imaginary options");
      sep = 0.0;
      break;
    }

  /* Select cases (dpe or floating point)
   * First normalize the polynomial (only the float version) */
  rdpe_div(tmp, max_coeff, min_coeff);
  rdpe_mul_eq_d(tmp, (double) (n + 1));
  rdpe_mul_eq(tmp, rdpe_mind);
  rdpe_div_eq(tmp, rdpe_maxd);
  if (rdpe_lt(tmp, rdpe_one)) {
    /* if  (n+1)*max_coeff/min_coeff < dhuge/dtiny -  float case */
    *which_case = 'f';
    rdpe_mul_eq(min_coeff, max_coeff);
    rdpe_mul(tmp, rdpe_mind, rdpe_maxd);
    rdpe_div(min_coeff, tmp, min_coeff);
    rdpe_sqrt_eq(min_coeff);
    /* min_coeff = sqrt(dhuge*dtiny/(min_coeff*max_coeff)) */
    for (i = 0; i <= n; i++) {
      rdpe_mul(tmp, dap[i], min_coeff);
      fap[i] = rdpe_get_d(tmp);
      if (data_type[1] == 'c') {
	cdpe_mul_e(ctmp, dpc[i], min_coeff);
	cdpe_get_x(fpc[i], ctmp);
      } else {
	/* fpr(i)=dpr(i)*min_coeff !! DARIO riattiva dopo impl. caso reale */
	cdpe_mul_e(ctmp, dpc[i], min_coeff);
	cdpe_get_x(fpc[i], ctmp);
      }
    }
  } else
    *which_case = 'd';
}

/*********************************************************
*      SUBROUTINE COMPUTE_SEP                            *
*********************************************************/
void
compute_sep(void)
{
  sep = n * lmax_coeff;
  sep = -sep - n * (1 + log((double) n)) / LOG2;
  if (DOLOG)
    fprintf(logstr, "Sep = %f\n", sep);
}
