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
 *                       IMPROVE
 ***********************************************************
 Improve all the approximations up to prec_out digits
 For each approximation compute the value of sigma such that
    e_j < e_0 * sigma **(2**j)  sigma=k/(k-1)=1/(1-t), k=1/t,
    t = Min_j |root(i)-root(j)|-rad(j)
 Then compute the number of digits needed for the j-th
 iteration i.e., digits=log(e_j/|x|)+cond, where
 log(e_j/|x|) = f+g**(2**j), ...
 cond = log(rad/eps), 
 cond ~ norm(p)(1+|root(i)|/(a_n prod_{j != i}|root(i)-root(j)|
 and cond ~ radius(i)/(eps |root(i)|) for user-defined polyn.
 mpwp denotes the number of bits of the current working
 precision.
 ************************************************************/
void
mps_improve(mps_status* st)
{
  int i, j, k, m;
  long mpnb_in, mpnb_out;
  tmpc_t mtmp;
  tmpc_t nwtcorr;
  cdpe_t ctmp;
  rdpe_t tmp, t, s, sigma, newrad, oldrad, abroot;
  double f, g, cnd;
  boolean again;

  if (st->DOLOG)
    fprintf(st->logstr, "Refining the roots ...\n");

  /* == 1 ==
   * compute the number mpnb_in of bits
   * corresponding to the given input precision. 
   * Set mpnb_in=0 if the input precision is infinite (prec_in=0) */
  if (st->prec_in == 0)
    mpnb_in = 0;
  else
    mpnb_in = (long) (st->prec_in * LOG2_10 + log(4.0 * st->n) / LOG2);
  mpnb_out = (long) (st->prec_out * LOG2_10);

  /* == 2  ==
   * compute the coefficients of the polynomial as mpc_t with mpnb_in bits
   * only if the polynomial is not assigned as a straight line program and
   * the input precision is not infinite. */
  if (mpnb_in != 0)
    mp_set_prec(s, mpnb_in);

  //  tmpc_init2(mtmp, mpwp); /* puo' essere settato a precisione minima */
  tmpc_init2(mtmp, mpnb_out*2); /* puo' essere settato a precisione minima */
  tmpc_init2(nwtcorr, mpnb_out*2);

  if (st->prec_in != 0 && st->data_type[0] != 'u')
    mps_prepare_data(st, mpnb_in);
  else {
    mp_set_prec(s, mpnb_out * 2);
    mps_prepare_data(st, mpnb_out * 2);
  }

  /* == 3 == 
   * scan the approximations to apply Newton's iterations */
  for (i = 0; i < st->n; i++) {
    if (st->DOLOG)
      fprintf(st->logstr, "root %d\n", i);
    if (st->status[i][0] != 'i' || st->status[i][2] == 'o')
      continue;			/* Do not refine approximated roots */

    /*  == 3.1 ==
     * for data_type[0]='d' compute  t=Min_j |root(i)-root(j)|-rad(j)-rad(i) 
     * otherwise set t=5*n*rad[i] since the root is Newton-isolated.
     * This allows us to remove an O(n^2) complexity  */

    if (st->data_type[0] == 's')
      rdpe_mul_d(t, st->drad[i], 5.0 * st->n);
    else {
      k = i + 1;
      if (i == st->n - 1)
	k = 0;
      mpc_sub(mtmp, st->mroot[k], st->mroot[i]);
      mpc_get_cdpe(ctmp, mtmp);
      cdpe_mod(t, ctmp);
      rdpe_sub_eq(t, st->drad[k]);
      rdpe_sub_eq(t, st->drad[i]);
      for (j = 0; j < st->n; j++)
	if (j != i) {
	  mpc_sub(mtmp, st->mroot[j], st->mroot[i]);
	  mpc_get_cdpe(ctmp, mtmp);
	  cdpe_mod(tmp, ctmp);
	  rdpe_sub_eq(tmp, st->drad[i]);
	  rdpe_sub_eq(tmp, st->drad[j]);
	  if (rdpe_gt(t, tmp))
	    rdpe_set(t, tmp);
	}
    }

    /*  == 3.2 ==
     * compute an  estimate of the condition number in terms of bits
     * as log_2(rad/(4*n*epsilon*|x|))       */

    rdpe_mul_d(tmp, st->drad[i], 4.0 * st->n);
    mpc_get_cdpe(ctmp, st->mroot[i]);
    cdpe_mod(abroot, ctmp);
    rdpe_div(tmp, tmp, abroot);
    cnd = st->rootwp[i] + rdpe_log(tmp) / LOG2 + 1;

    /* then evaluate the number of bits g,f */
    rdpe_div(t, st->drad[i], t);
    rdpe_mul_eq_d(t, (double) st->n - 1);
    rdpe_sub(s, rdpe_one, t);
    rdpe_div(sigma, t, s);
    g = -rdpe_log(sigma) / LOG2;
    rdpe_set(tmp, abroot);
    rdpe_mul_eq(tmp, sigma);
    rdpe_div(tmp, st->drad[i], tmp);
    f = -rdpe_log(tmp) / LOG2;

    /* evaluate the upper bound m to the number of iterations
     * needed to reach the desired precision */
    m = (int) (log((mpnb_out - f) / g) / LOG2) + 1;

    /*  == 4 ==      Start Newton */

    rdpe_set(oldrad, st->drad[i]);
    for (j = 1; j <= m; j++) {
      if (st->DOLOG)
	fprintf(st->logstr, "iter= %d\n", j);
      g *= 2;
      st->mpwp = (long) (f + g + cnd);
      if (st->mpwp >= mpnb_in && mpnb_in != 0)
	st->mpwp = mpnb_in;

      tmpc_clear(nwtcorr);
      tmpc_init2(nwtcorr, st->mpwp);

      mp_set_prec(s, st->mpwp);
      mps_prepare_data(st, st->mpwp);
      if (st->data_type[0] != 'u')
	mps_mnewton(st, st->n, st->mroot[i], st->drad[i], 
		    nwtcorr, st->mfpc, st->mfppc, st->dap, st->spar, &again);
      else
	mps_mnewton_usr(st, st->mroot[i], st->drad[i], nwtcorr, &again);
      mpc_sub_eq(st->mroot[i], nwtcorr);

      /* correct radius, since the computed one is referred to the previous
       * approximation. Due to the quadratic convergence the new approximation
       * the radius is bounded by 2^(-g-f+1)
       */
      rdpe_set_2dl(newrad, 1.0, (long) (-g - f + 1));
      rdpe_set(tmp, abroot);
      rdpe_mul_eq(newrad, tmp);
      rdpe_mul_eq(tmp, st->eps_out);

      if (rdpe_eq(st->drad[i], rdpe_zero))
	rdpe_set(st->drad[i], newrad);
      if (rdpe_lt(newrad, st->drad[i]))
	rdpe_set(st->drad[i], newrad);

      if (rdpe_lt(st->drad[i], tmp) || st->mpwp == mpnb_in)
	break;			/* loop1 */
    }

    /* update the record working precision for root i */
    st->rootwp[i] = mpc_get_prec(st->mroot[i]);
  }

  tmpc_clear(nwtcorr);
  tmpc_clear(mtmp);
}

