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

#include <mps/core.h>
#include <mps/rursolve.h>

/*********************************************************
*                   PROCEDURE HORNER                     *
**********************************************************/
void
mps_horner (mps_status * s, mpc_t res, int *dprec, int *iprec, int deg, int i)
{
  int j, k;
  tmpc_t x, y;
  cdpe_t ctmp;
  rdpe_t tmp, tmp1, tmp2;

  /* compute the precision of mroot[i] in bits */
  mpc_get_cdpe (ctmp, s->mroot[i]);
  cdpe_mod (tmp, ctmp);
  rdpe_div (tmp1, s->drad[i], tmp);
  k = -rdpe_Esp (tmp1);
  *iprec = k;

  tmpc_init2 (x, k + 4);
  tmpc_init2 (y, k + 4);

  mpf_set_z (mpc_Re (y), mpdemo[deg]);
  mpf_set_ui (mpc_Im (y), 0);
  for (j = deg - 1; j >= 0; j--)
    {
      mpc_mul_eq (y, s->mroot[i]);
      mpf_set_z (mpc_Re (x), mpdemo[j]);
      mpf_set_ui (mpc_Im (x), 0);
      mpc_add_eq (y, x);
    }

  rdpe_set (tmp2, s->dap1[deg]);
  for (j = deg - 1; j >= 0; j--)
    {
      rdpe_mul_eq (tmp2, tmp);
      rdpe_add_eq (tmp2, s->dap1[j]);
    }

  rdpe_mul_eq_d (tmp1, 4.0 * deg);
  rdpe_mul_eq (tmp2, tmp1);
  mpc_get_cdpe (ctmp, y);
  cdpe_mod (tmp, ctmp);
  rdpe_div_eq (tmp2, tmp);
  *dprec = -rdpe_Esp (tmp2);

  if (deg == 0)
    *dprec = *iprec;
  mpc_set (res, y);

  tmpc_clear (y);
  tmpc_clear (x);
}

/*********************************************************
*                   PROCEDURE REFINE                     *
**********************************************************/
void
mps_refine (mps_status * s, int i, long prec)
{
  int j, k, m;
  tmpc_t mtmp;
  mpc_t nwtcorr;
  cdpe_t ctmp;
  rdpe_t tmp, t, cond, norm, u, sigma;
  double f, g, cnd;
  mps_boolean aga;

  tmpc_init2 (mtmp, prec);
  mpc_init2 (nwtcorr, prec);
  rdpe_set_2dl (s->eps_out, 1, -prec);

  /* ======
   * compute the infinite norm of the polynomial p/p_{n+1}
   * ======= */
  rdpe_set (norm, s->dap[s->n]);
  for (j = 0; j < s->n; j++)
    if (rdpe_lt (norm, s->dap[j]))
      rdpe_set (norm, s->dap[j]);
  rdpe_div_eq (norm, s->dap[s->n]);

  /* Do not refine nonisolated roots */
  if (s->status[i][0] != 'i')
    {
      fprintf (s->logstr, "Warning nonisolated root for i=%d\n", i);
      return;
    }

  /*  =====================================
   * Compute the number of iterations needed by Newton 
   * ===================================== */

  /* first:
   * evaluate t=Min_j |root(i)-root(j)|-rad(j)
   */

  k = i + 1;
  if (i == s->n - 1)
    k = 0;

  mpc_sub (mtmp, s->mroot[k], s->mroot[i]);
  mpc_get_cdpe (ctmp, mtmp);
  cdpe_mod (t, ctmp);
  rdpe_sub_eq (t, s->drad[k]);
  rdpe_set (cond, rdpe_one);
  for (j = 0; j < s->n; j++)
    if (j != i)
      {
	mpc_sub (mtmp, s->mroot[j], s->mroot[i]);
	mpc_get_cdpe (ctmp, mtmp);
	cdpe_mod (tmp, ctmp);
	rdpe_sub_eq (tmp, s->drad[i]);
	if (rdpe_gt (t, tmp))
	  rdpe_set (t, tmp);
	rdpe_sub_eq (tmp, s->drad[j]);
	rdpe_mul_eq (cond, tmp);
      }

  /* second:
   * compute an  estimate of the condition number as 
   * cond= (1+|x|**n)*norm(p)/prod{i != j}|root(i)-root(j)|
   */

  mpc_get_cdpe (ctmp, s->mroot[i]);
  cdpe_mod (tmp, ctmp);
  rdpe_mul_eq (cond, tmp);
  rdpe_pow_si (u, tmp, s->n);
  rdpe_mul_eq (u, norm);
  rdpe_add_eq (u, norm);
  rdpe_div (cond, u, cond);
  cnd = rdpe_log (cond) / LOG2 + 1;

  /* third:  
   * evaluate the number of bits g,f
   */

  rdpe_div (t, s->drad[i], t);
  rdpe_mul_eq_d (t, (double) s->n - 1);
  rdpe_sub (u, rdpe_one, t);
  rdpe_div (sigma, t, u);
  g = -rdpe_log (sigma) / LOG2;
  mpc_get_cdpe (ctmp, s->mroot[i]);
  cdpe_mod (tmp, ctmp);
  rdpe_mul_eq (tmp, sigma);
  rdpe_div (tmp, s->drad[i], tmp);
  f = -rdpe_log (tmp) / LOG2;

  /* finally
   * evaluate the upper bound m to the number of iterations
   * needed to reach the desired precision
   */

  m = (int) (log ((prec - f) / g) / LOG2) + 1;

  /* ==========================
   * Start Newton    
   * ========================== */
  for (j = 1; j <= m; j++)
    {
      if (s->DOLOG)
	fprintf (s->logstr, "iter= %d\n", j);
      g *= 2;
      s->mpwp = (long) (f + g + cnd);
      mpc_set_prec (nwtcorr, s->mpwp);
      mps_mp_set_prec (s, s->mpwp);
      mps_prepare_data (s, s->mpwp);
      mps_mnewton (s, s->n, s->mroot[i], s->drad[i], nwtcorr, s->mfpc,
		   s->mfppc, s->dap, s->spar, &aga, 0);
      mpc_sub_eq (s->mroot[i], nwtcorr);

      /* correct radius, since the computed one is referred to the previous
       * approximation. Due to the quadratic convergence the new approximation
       * has at least g/2 more correct bits. */
      rdpe_set_2dl (tmp, 1.0, (long int) (-g / 2.0));
      rdpe_mul_eq (s->drad[i], tmp);
      mpc_get_cdpe (ctmp, s->mroot[i]);
      cdpe_mod (tmp, ctmp);
      rdpe_mul_eq (tmp, s->eps_out);
      if (rdpe_lt (s->drad[i], tmp))
	break;			/* loop1 */
    }

  /* Check for zero radius */
  if (rdpe_eq (s->drad[i], rdpe_zero))
    {
      mpc_get_cdpe (ctmp, s->mroot[i]);
      cdpe_mod (tmp, ctmp);
      rdpe_mul_eq (tmp, s->eps_out);
      rdpe_set (s->drad[i], tmp);
    }

  mpc_clear (nwtcorr);
  tmpc_clear (mtmp);
}
