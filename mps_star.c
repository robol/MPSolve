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

/**
 * @file 
 * @brief Routines to compute starting approximations
 * for the algorithm
 *
 */

#include "mps.h"

static const double pi2 = 6.283184;

/* forward declaration */
void raisetemp(unsigned long int digits);
void raisetemp_raw(unsigned long int digits);

/**
 * @brief Compute new starting approximations to the roots
 * of the polynomial \f$p(x)\f$ having coefficients of modulus apoly.
 *
 * Computations is done by
 * means of the Rouche'-based criterion of Bini (Numer. Algo. 1996). 
 * The program can compute all the approximations
 * (if \f$n\f$ is the degree of \f$p(x)\f$) or it may compute the
 * approximations of the cluster of index <code>i_clust</code>
 * The status vector is changed into <code>'o'</code> for the components
 * that belong to a cluster with relative radius less than <code>eps</code>.
 * The status vector is changed into <code>'x'</code> for the components that
 * cannot be represented as double.
 *
 * @param n number of roots in the cluster.
 * @param i_clust index of cluster to analyze.
 * @param clust_rad radius of cluster.
 * @param g gravity center of the cluster.
 * @param eps a double that represent the maximum value
 * of relative radius (with respect to <code>g</code>) of
 * roots whose status must be set to <code>o</code>.
 * @param fap array of moduli of the coefficients as double.
 *
 * @see status
 */
void
fstart(int n, int i_clust, double clust_rad, double g, rdpe_t eps,
       double fap[])
{
  const double  big = DBL_MAX,   small = DBL_MIN;
  const double xbig = log(big), xsmall = log(small);

  int i, iold, j, jj, l, ni, nzeros;
  double sigma, th, ang, temp, r;
  rdpe_t tmp;

  if (random_seed)
    sigma = drand();
  else
    sigma=0.1;

  ni = 0;
  nzeros = 0;
  r = 0.0;

  /* In the case of user-defined polynomial choose as starting
   * approximations equispaced points in the unit circle.  */
  if (data_type[0] == 'u') {
    ang = pi2 / n;
    for (i = 0; i < n; i++)
      cplx_set_d(froot[i], cos(ang * i + sigma), sin(ang * i + sigma));
    return;
  }

  /* In the general case apply the Rouche-based criterion */

  /********************************************** 
    check for possible null entries in the trailing
    coefficients only in the case where the polynomial
    has been shifted in g, replace null coefficients
    with small numbers according to the working precision
    and to the number of null coefficients
    **********************************************/

  if (g != 0.0) {
    for (i = 0; i <= n; i++)
      if (fap[i] != 0.0) {
	ni = i;
	break;
      }
    if (ni == 0)
      temp = 2 * xsmall;
    else
      temp = log(fap[ni]) + ni * (log(DBL_EPSILON) + log(g * ni * 10.0));
  } else
    temp = 2 * xsmall;

  for (i = 0; i <= n; i++)
    if (fap[i] != 0.0)
      fap2[i] = log(fap[i]);
    else
      fap2[i] = temp;

  /* compute the convex hull */
  fconvex(n, fap2);

  /* compute the radii of the circles containing starting approximations  */
  iold = 0;
  th = pi2 / n;
  for (i = 1; i <= n; i++)
    if (h[i]) {
      nzeros = i - iold;
      temp = (fap2[iold] - fap2[i]) / nzeros;
      /* if the radius is too small to be represented as double, set it
       * to the minimum  representable double */
      if (temp < xsmall)	/* if (temp < MAX(xsmall, -xbig)) DARIO Giugno 23 */
	r = DBL_MIN;		/* r = small; */

      /* if the radius is too big to be represented as double, set it
       * to the maximum representable double */
      if (temp > xbig)
	r = DBL_MAX;		/* big;   DARIO Giugno 23 */

      /* if the radius is representable as double, compute it    */
      if ((temp <= xbig) && (temp > xsmall))
	/* if ((temp <= xbig) && (temp > MAX(-xbig, xsmall))) DARIO Giugno 23 */
	r = exp(temp);

      /* if the radius is greater than the radius of the cluster
       * set the radius equal to the radius of the cluster */
      if (clust_rad != 0 && r > clust_rad)
	r = clust_rad;

      /* Choose starting values for root */
      ang = pi2 / nzeros;
      for (j = iold; j <= i - 1; j++) {
	if (g != 0.0)
	  l = clust[punt[i_clust] + j];
	else
	  l = j;
	jj = j - iold;

	/* if the radius reaches extreme values then set the components
	 * of status, corresponding to approximation which fall out the 
	 * representable range, to 'x' (out)    */
	if ((r == DBL_MIN) || (r == DBL_MAX))
	  /* if ((r == small) || (r == big)) DARIO Giugno 23 */
	  status[l][0] = 'x';
	cplx_set_d(froot[l], r * cos(ang * jj + th * i + sigma),
		   r * sin(ang * jj + th * i + sigma));
      }
      iold = i;
    }
  /* If the new radius of the cluster is relatively smaller, then
   * set the status component equal to 'o' (output) */
  if (g != 0.0) {
    rdpe_mul_d(tmp, eps, g);
    if (r * nzeros <= rdpe_get_d(tmp))
      for (j = 0; j < punt[i_clust + 1] - punt[i_clust]; j++) {
	l = clust[punt[i_clust] + j];
	status[l][0] = 'o';
	frad[l] = r * nzeros;
      }
  }
}

/*********************************************************
*              SUBROUTINE DSTART                         *
**********************************************************
 Compute new starting approximations to the roots of the
 polynomial p(x) having coefficients of modulus apoly, by
 means of the Rouche'-based criterion of Bini (Numer. Algo. 1996). 
 The program can compute all the approximations
 (if n is the degree of p(x)) or it may compute the
 approximations of the cluster of index i_clust
 The status vector is changed into 'o' for the components
 that belong to a cluster with relative radius less than eps.
 The status vector is changed into 'f' for the components
 that cannot be represented as dpe.
 ***********************************************/
void
dstart(int n, int i_clust, rdpe_t clust_rad,
       rdpe_t g, rdpe_t eps, rdpe_t dap[])
{
  int l, i, j, jj, iold, ni = 0, nzeros = 0;
  rdpe_t r, tmp, tmp1;
  double sigma, th, ang, xbig, xsmall, temp;
  boolean flag = false;

  if (random_seed) 
    sigma = drand();
  else
    sigma=0.1;

  /* In the case of user-defined polynomial choose as starting
   * approximations equispaced points in the unit circle. */
  if (data_type[0] == 'u') {
    ang = pi2 / n;
    for (i = 0; i < n; i++)
      cdpe_set_d(droot[i], cos(ang * i + sigma), sin(ang * i + sigma));
    return;
  }
  
  /* In the general case apply the Rouche-based criterion */

  xsmall = rdpe_log(RDPE_MIN);
  xbig = rdpe_log(RDPE_MAX);

  /* check if it is the case dpe_after_float, in this case set flag=true  */
  for (i = 0; i < n; i++) {
    flag = (status[i][0] == 'x');
    if (flag)
      break;
  }

  /* check for possible null entries in the trailing coefficients
   * only in the case where the polynomial has been shifted in g
   * replace null coefficients with small numbers according to
   * the working precision and to the number of null coefficients */

  if (rdpe_ne(g, rdpe_zero)) {
    for (i = 0; i <= n; i++)
      if (rdpe_ne(dap[i], rdpe_zero)) {
	ni = i;
	break;
      }
    if (ni == 0)
      temp = -2.0 * (LONG_MAX * LOG2);
    else {
      /* temp = log(dap[ni])+ni*(log(DBL_EPSILON)+log(g*ni*10.d0)) */
      temp = ni * 10.0;
      rdpe_mul_d(tmp, g, temp);
      temp = rdpe_log(tmp);
      temp += log(DBL_EPSILON);
      temp *= ni;
      temp += rdpe_log(dap[ni]);
    }
  } else
    temp = -2.0 * (LONG_MAX * LOG2);
  for (i = 0; i <= n; i++)
    if (rdpe_ne(dap[i], rdpe_zero))
      fap2[i] = rdpe_log(dap[i]);
    else
      fap2[i] = temp;

  /* compute the convex hull */
  fconvex(n, fap2);
  th = pi2 / n;

  /* Scan all the vertices of the convex hull   */
  iold = 0;
  for (i = 1; i <= n; i++)
    if (h[i]) {
      nzeros = i - iold;
      temp = (fap2[iold] - fap2[i]) / nzeros;

      /* if the radius is too small or too big to be represented as dpe, 
       * output a warning message */
      if (temp < xsmall) {
	rdpe_set(r, RDPE_MIN);
	if (DOLOG) {
	  fprintf(logstr, "Warning: Some zeros are too small to be\n");
	  fprintf(logstr, " represented as cdpe, they are replaced by\n");
	  fprintf(logstr, " small numbers and the status is set to 'F'.\n");
	}
      }
      if (temp > xbig) {
	rdpe_set(r, RDPE_MAX);
	if (DOLOG) {
	  fprintf(logstr, "Warning: Some zeros are too big to be\n");
	  fprintf(logstr, " represented as cdpe, they are replaced by\n");
	  fprintf(logstr, " big numbers and the status is set to 'F'.\n");
	}
      }

      /* if the radius is representable as dpe, compute it */
      if ((temp <= xbig) && (temp >= xsmall)) {
	rdpe_set_d(r, temp);
	rdpe_exp_eq(r);
      }

      /* if the radius is greater than the radius of the cluster
       * set the radius equal to the radius of the cluster */
      if (rdpe_ne(g, rdpe_zero) && rdpe_gt(r, clust_rad))
	rdpe_set(r, clust_rad);

      ang = pi2 / nzeros;
      for (j = iold; j < i; j++) {
	if (rdpe_ne(g, rdpe_zero))
	  l = clust[punt[i_clust] + j];
	else
	  l = j;

	jj = j - iold;

	/* If dpe_after_float (i.e., flag is true) recompute the starting
	 * values of only the approximations falling out of the range */
	if (flag) {
	  if (status[l][0] == 'x') {
	    cdpe_set_d(droot[l], cos(ang * jj + th * i + sigma),
		       sin(ang * jj + th * i + sigma));
	    cdpe_mul_eq_e(droot[l], r);
	    status[l][0] = 'c';
	    /*#G 27/4/98 if (rdpe_eq(r, big) || rdpe_eq(r, small)) */
	    if (rdpe_eq(r, RDPE_MIN) || rdpe_eq(r, RDPE_MAX))
	      status[l][0] = 'f';
	  }
	} else {
	  /* else compute all the initial approximations */
	  cdpe_set_d(droot[l], cos(ang * jj + th * i + sigma),
		     sin(ang * jj + th * i + sigma));
	  cdpe_mul_eq_e(droot[l], r);
	  /*#G 27/4/98 if (rdpe_eq(r, big) || rdpe_eq(r, small)) */
	  if (rdpe_eq(r, RDPE_MIN) || rdpe_eq(r, RDPE_MAX))
	    status[l][0] = 'f';
	}
      }
      iold = i;
    }
  /* If the new radius of the cluster is relatively small, then
   * set the status component equal to 'o' (output) */
  if (rdpe_ne(g, rdpe_zero)) {
    rdpe_mul(tmp, g, eps);
    rdpe_mul_d(tmp1, r, (double) nzeros);
    if (rdpe_lt(tmp1, tmp))
      for (j = 0; j <= punt[i_clust + 1] - punt[i_clust]; j++) {
	l = clust[punt[i_clust] + j];
	status[l][0] = 'o';
	rdpe_set(drad[l], tmp1);
      }
  }
}

/***********************************************************
*            SUBROUTINE MSTART                             *
***********************************************************/
void
mstart(int n, int i_clust, rdpe_t clust_rad, rdpe_t g,
       rdpe_t dap[])
{
  int i, j, jj, iold, l, nzeros;
  double sigma, ang, th, xbig, xsmall, temp;
  rdpe_t r, big, small, rtmp1, rtmp2;
  cdpe_t ctmp;

  if (random_seed) 
    sigma = drand();
  else
    sigma=0.1;

  nzeros = 0;
  temp = 0.0;

  /* In the general case apply the Rouche-based criterion */

  xsmall = rdpe_log(RDPE_MIN);
  xbig = rdpe_log(RDPE_MAX);
  rdpe_set(small, RDPE_MIN);
  rdpe_set(big, RDPE_MAX);

  if (rdpe_eq(dap[0], rdpe_zero))
    fap2[0] = -mpwp * LOG2;

  /*  check for possible null entries in the trailing coefficients */
  for (i = 0; i <= n; i++)
    if (rdpe_ne(dap[i], rdpe_zero))
      fap2[i] = rdpe_log(dap[i]);
    else
      fap2[i] = fap2[0];

  /* compute the convex hull */
  fconvex(n, fap2);
  th = pi2 / n;

  /* Scan all the vertices of the convex hull */
  iold = 0;
  for (i = 1; i <= n; i++)
    if (h[i]) {
      nzeros = i - iold;
      temp = (fap2[iold] - fap2[i]) / nzeros;
      /* if the radius is too small or too big to be represented as dpe, 
       * output a warning message */
      if (temp < xsmall) {
	rdpe_set(r, small);
	if (DOLOG) {
	  fprintf(logstr, "Warning: Some zeros are too small to be\n");
	  fprintf(logstr, " represented as cdpe, they are replaced by\n");
	  fprintf(logstr, " small numbers and the status is set to 'F'.\n");
	}
      }
      if (temp > xbig) {
	rdpe_set(r, big);
	if (DOLOG) {
	  fprintf(logstr, "Warning: Some zeros are too big to be\n");
	  fprintf(logstr, " represented as cdpe, they are replaced by\n");
	  fprintf(logstr, " big numbers and the status is set to 'F'.\n");
	}
      }
      /* if the radius is representable as dpe, compute it */

      if (temp <= xbig && temp >= xsmall) {
	rdpe_set_d(r, temp);
	rdpe_exp_eq(r);
      }
      /* if the radius is greater than the radius of the cluster
       * set the radius equal to the radius of the cluster */
      if (rdpe_gt(r, clust_rad)) 
	rdpe_set(r, clust_rad); 
      ang = pi2 / nzeros;

      /* Compute the initial approximations */
      for (j = iold; j < i; j++) {
	jj = j - iold;
	l = clust[punt[i_clust] + j];
	cdpe_set_d(ctmp, cos(ang * jj + th * i + sigma),
		   sin(ang * jj + th * i + sigma));
	cdpe_mul_eq_e(ctmp, r);
	cdpe_set(droot[l], ctmp);
	if (rdpe_eq(r, big) || rdpe_eq(r, small))
	  status[l][0] = 'f';
      }
      iold = i;
    }
  /* If the new radius of the cluster is relatively small, then
   * set the status component equal to 'o' (output)
   * and set the corresponding radius */

  rdpe_set(rtmp1, r);
  rdpe_mul_eq_d(rtmp1, (double) nzeros);
  rdpe_set(rtmp2, g);
  rdpe_mul_eq(rtmp2, eps_out);
  if (rdpe_le(rtmp1, rtmp2))
    for (j = 0; j < punt[i_clust + 1] - punt[i_clust]; j++) {
      l = clust[punt[i_clust] + j];
      status[l][0] = 'o';
      rdpe_mul_d(drad[l], r, (double) nzeros);
    }
  rdpe_set(clust_rad,r); 
}

/***********************************************************
*                     Subroutine FRESTART                  *
************************************************************
 The program scans the existing clusters and  selects the
 ones where shift in the gravity center must be done.
 Then computes the gravity center g, performs the shift of
 the variable and compute new starting approximations in the
 cluster.
 The components of the vector status(:,1) are set to 'c'
 (i.e., Aberth's iteration must be applied) if the cluster
 intersects the origin (in this case shift is not applied),
 or if new starting approximations have been selected.
 The gravity center g is choosen as a zero of the (m-1)-st
 derivative of the polynomial in the cluster, where m=m_clust
 is the multiplicity of the cluster. 

 Shift in g is perfomed if
        status=(c*u) for goal=count
        status=(c*u) or status=(c*i) for goal=isolate/approximate
 To compute g, first compute the weighted mean (super center sc)
 of the approximations in the cluster, where the weight are the
 radii, then compute the radius (super radius sr) of the disk
 centered in the super center containing all the disks of the cluster.
 Apply few steps of Newton's iteration to the (m-1)-st derivative
 of the polynomial starting from the super center and obtain the
 point g where to shift the variable.
 If g is outside the super disk of center sc and radius sr
 output a warning message.
 ***********************************************************/
void
frestart(void)
{
  int i, k, j, l, jj;
  double sr, sum, rad, rtmp, rtmp1;
  cplx_t sc, g, corr, ctmp;
  boolean tst, cont;

  /* For user's polynomials skip the restart stage (not yet implemented) */
  if (data_type[0] == 'u')
    return;

  /* scan the existing clusters and  select the ones where shift in
   * the gravity center must be done. tst=true means do not perform shift */
  for (i = 0; i < nclust; i++) {	/* loop1: */
    if ((punt[i + 1] - punt[i]) == 1)
      continue;
    tst = true;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* looptst : */
      l = clust[punt[i] + j];
      if (!again[l])
	goto loop1;
      if (goal[0] == 'c') {
	if (status[l][0] == 'c' && status[l][2] == 'u') {
	  tst = false;
	  break;
	}
      } else if ((status[l][0] == 'c' && status[l][2] == 'u')
		 || (status[l][0] == 'c' && status[l][2] == 'i')) {
	tst = false;
	break;
      }
    }				/* for */
    if (tst)
      goto loop1;

    /* Compute super center sc and super radius sr */
    sum = 0.0;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      sum += frad[l];
    }
    cplx_set(sc, cplx_zero);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      cplx_mul_d(ctmp, froot[l], frad[l]);
      cplx_add_eq(sc, ctmp);
    }
    cplx_div_eq_d(sc, sum);

    sr = 0.0;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      cplx_sub(ctmp, sc, froot[l]);
      sr = MAX(sr, frad[l] + cplx_mod(ctmp));
    }

    /* Check the relative width of the cluster
     * If it is greater than 1 do not shift
     * and set status(:1)='c' that means
     * keep iterating Aberth's step. */

    if (sr > cplx_mod(sc)) {
      for (j = punt[i]; j < punt[i + 1]; j++)
	status[clust[j]][0] = 'c';
      if (DOLOG)
	fprintf(logstr, "     FRESTART: cluster rel. large: skip to the next component\n");
      goto loop1;
    }
    
    /* Now check the Newton isolation of the cluster */

    for (k = 0; k < nclust; k++)
      if (k != i) {
	for (j = 0; j < punt[k + 1] - punt[k]; j++) {
	  cplx_sub(ctmp, sc, froot[clust[punt[k] + j]]);
	  rtmp = cplx_mod(ctmp);
	  rtmp1 = (sr + frad[clust[punt[k] + j]]) * 5 * n;
	  if (rtmp < rtmp1) {
	    for (jj = punt[i]; jj < punt[i + 1]; jj++)
	      status[clust[jj]][0] = 'c';
	    if (DOLOG) {
	      fprintf(logstr, "Cluster not Newton isolated:");
	      fprintf(logstr, "  skip to the next component\n");
	    }
	    goto loop1;
	  }
	}
      }
    /* Compute the coefficients of the derivative of p(x) having order
     * equal to the multiplicity of the cluster -1. */
    sum = 0.0;
    for (j = 0; j <= n; j++) {
      sum += cplx_mod(fpc[j]);
      cplx_set(fppc[j], fpc[j]);
    }
    for (j = 1; j < punt[i + 1] - punt[i]; j++) {
      for (k = 0; k <= n - j; k++)
	cplx_mul_d(fppc[k], fppc[k + 1], (double) (k + 1));
    }
    for (j = 0; j < n - (punt[i + 1] - punt[i]) + 2; j++)
      fap1[j] = cplx_mod(fppc[j]);

    /* Apply at most max_newt_it steps of Newton's iterations
     * to the above derivative starting from the super center
     * of the cluster. */
     
    cplx_set(g, sc);
    for (j = 0; j < max_newt_it; j++) {		/* loop_newt: */
      rad = 0.0;
      fnewton(n - (punt[i + 1] - punt[i]) + 1, g,
	      &rad, corr, fppc, fap1, &cont);
      cplx_sub_eq(g, corr);
      if (!cont)
	break;
    }
    if (j == max_newt_it) {
      if (DOLOG)
	fprintf(logstr, "Exceeded maximum Newton iterations in frestart\n");
      return;
    }
    cplx_sub(ctmp, sc, g);
    if (cplx_mod(ctmp) > sr) {
      if (DOLOG)
	fprintf(logstr, "The gravity center falls outside the cluster\n");
      return;
    }
    /* Compute the coefficients of the shifted polynomial p(x+g)
     * and compute new starting approximations
     * First check if shift may cause overflow, in this case skip
     * the shift stage */

    if (n * log(cplx_mod(g)) + log(sum) > log(DBL_MAX))
      goto loop1;
    if (DOLOG)
      fprintf(logstr, "      FRESTART:  fshift\n");
    fshift(punt[i + 1] - punt[i], i, sr, g, eps_out);
    rtmp = cplx_mod(g);
    rtmp *= DBL_EPSILON * 2;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      /* Choose as new incl. radius 2*multiplicity*(radius of the circle) */
      frad[l] = 2 * (punt[i + 1] - punt[i]) * cplx_mod(froot[l]);
      cplx_add_eq(froot[l], g);
      if (frad[l] < rtmp)	/* DARIO* aggiunto 1/5/97 */
	frad[l] = rtmp;
    }
  loop1:;
  }
}

/*************************************************************
*                     SUBROUTINE DRESTART                    *
**************************************************************
 The program scans the existing clusters and  selects the ones 
 where shift in the gravity center must be done.
 Then computes the gravity center g, performs the shift of the variable
 and compute new starting approximations in the cluster.
 The components of the vector err are set to true (i.e., Aberth's
 iteration must be applied) if the cluster intersects the origin
 (in this case shift is not applied),
 or if new starting approximations have been selected.
 The gravity center g is choosen as a zero of the (m-1)-st derivative
 of the polynomial in the cluster, where m=mclust is the multiplicity
 of the cluster. 

 Shift in g is perfomed if
        status=(c*u) for goal=count
        status=(c*u) or status=(c*i) for goal=isolate/approximate
 To compute g, first compute the weighted mean (super center sc)
 of the approximations in the cluster, where the weight are the radii,
 then compute the radius (super radius sr) of the disk centered in the
 super center containing all the disks of the cluster.
 Apply few steps of Newton's iteration to the (m-1)-st derivative
 of the polynomial starting from the super center and obtain the point
 g where to shift the variable.
 If g is outside the super disk of center sc and radius sr output a warning
 message.
 ******************************************************************/
void
drestart(void)
{
  int i, k, j, l, jj;
  rdpe_t sr, sum, rad, rtmp, rtmp1;
  cdpe_t sc, g, corr, ctmp;
  boolean tst, cont;

  /*  For user's polynomials skip the restart stage (not yet implemented) */
  if (data_type[0] == 'u')
    return;

  for (i = 0; i < nclust; i++) {	/* loop1: */
    if ((punt[i + 1] - punt[i]) == 1)
      continue;
    tst = true;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* looptst: */
      l = clust[punt[i] + j];
      if (!again[l])
	goto loop1;
      if (goal[0] == 'c') {
	if (status[l][0] == 'c' && status[l][2] == 'u') {
	  tst = false;
	  break;
	}
      } else if ((status[l][0] == 'c' && status[l][2] == 'u')
		 || (status[l][0] == 'c' && status[l][2] == 'i')) {
	tst = false;
	break;
      }
    }				/* for */
    if (tst)
      goto loop1;

    /* Compute super center sc and super radius sr */
    rdpe_set(sum, rdpe_zero);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      rdpe_add_eq(sum, drad[l]);
    }
    cdpe_set(sc, cdpe_zero);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      cdpe_mul_e(ctmp, droot[l], drad[l]);
      cdpe_add_eq(sc, ctmp);
    }
    cdpe_div_eq_e(sc, sum);
    rdpe_set(sr, rdpe_zero);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      cdpe_sub(ctmp, sc, droot[l]);
      cdpe_mod(rtmp, ctmp);
      rdpe_add_eq(rtmp, drad[l]);
      if (rdpe_lt(sr, rtmp))
	rdpe_set(sr, rtmp);
    }

    /* Check the relative width of the cluster
     * If it is greater than 1 do not shift
     * and set statu(:1)='c' that means
     * keep iterating Aberth's step. */
    cdpe_mod(rtmp, sc);
    if (rdpe_gt(sr, rtmp)) {
      for (j = punt[i]; j < punt[i + 1]; j++) {
	status[clust[j]][0] = 'c';
	/* err(clust[j])=true  */
      }
      if (DOLOG)
	fprintf(logstr, "     DRESTART: cluster rel. large: skip to the next component\n");
      goto loop1;
    }
    /* Now check the Newton isolation of the cluster */

    for (k = 0; k < nclust; k++)
      if (k != i) {
	for (j = 0; j < punt[k + 1] - punt[k]; j++) {
	  cdpe_sub(ctmp, sc, droot[clust[punt[k] + j]]);
	  cdpe_mod(rtmp, ctmp);
	  rdpe_add(rtmp1, sr, drad[clust[punt[k] + j]]);
	  rdpe_mul_eq_d(rtmp1, 2.0 * n);     
	  if (rdpe_lt(rtmp, rtmp1)) {
	    for (jj = punt[i]; jj < punt[i + 1]; jj++)
	      status[clust[jj]][0] = 'c';
	    if (DOLOG) {
	      fprintf(logstr, "cluster not Newton isolated:");
	      fprintf(logstr, " skip to the next component\n");
	    }
	    goto loop1;
	  }
	}
      }
      
    /* Compute the coefficients of the derivative of p(x) having order
     * equal to the multiplicity of the cluster -1. */
     
    for (j = 0; j <= n; j++)
      cdpe_set(dpc2[j], dpc[j]);
    for (j = 1; j < punt[i + 1] - punt[i]; j++) {
      for (k = 0; k <= n - j; k++)
	cdpe_mul_d(dpc2[k], dpc2[k + 1], (double) (k + 1));
    }
    for (j = 0; j < n - (punt[i + 1] - punt[i]) + 2; j++)
      cdpe_mod(dap1[j], dpc2[j]);

    /* Apply at most max_newt_it steps of Newton's iterations
     * to the above derivative starting from the super center
     * of the cluster. */
     
    cdpe_set(g, sc);
    for (j = 0; j < max_newt_it; j++) {		/* loop_newt: */
      rdpe_set(rad, rdpe_zero);
      dnewton(n - (punt[i + 1] - punt[i]) + 1, g, rad,
	      corr, dpc2, dap1, &cont);
      cdpe_sub_eq(g, corr);
      if (!cont)
	break;
    }
    if (j == max_newt_it) {
      if (DOLOG)
	fprintf(logstr, "Exceeded maximum Newton iterations in frestart\n");
      return;
    }
    cdpe_sub(ctmp, sc, g);
    cdpe_mod(rtmp, ctmp);
    if (rdpe_gt(rtmp, sr)) {
      if (DOLOG)
	fprintf(logstr, "The gravity center falls outside the cluster\n");
      return;
    }
    /* Shift the variable and compute new approximations */
    if (DOLOG)
      fprintf(logstr, "      DRESTART:  dshift");
    dshift(punt[i + 1] - punt[i], i, sr, g, eps_out);
    cdpe_mod(rtmp, g);
    rdpe_mul_eq_d(rtmp, DBL_EPSILON * 2);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];

      /* Choose as new incl. radius 2*multiplicity*(radius of the circle) */
      cdpe_mod(drad[l], droot[l]);
      rdpe_mul_eq_d(drad[l], (double) (2 * (punt[i + 1] - punt[i])));
      cdpe_add_eq(droot[l], g);
      if (rdpe_lt(drad[l], rtmp))
	rdpe_set(drad[l], rtmp);
    }
  loop1:;
  }
}

/*************************************************************
*                     SUBROUTINE MRESTART                    *
*************************************************************/
void
mrestart(void)
{
  boolean tst, cont;
  int i, j, k, l, jj;
  rdpe_t sr, rad, rtmp, rtmp1, rtmp2;
  cdpe_t tmp;
  tmpf_t rea, srmp;
  tmpc_t sc, corr, temp;
  mpc_t g;
  
  /* For user's polynomials skip the restart stage (not yet implemented) */
  if (data_type[0] == 'u')
    return;

  tmpf_init2(rea, mpwp);
  tmpf_init2(srmp, mpwp);
  tmpc_init2(sc, mpwp);
  tmpc_init2(corr, mpwp);
  tmpc_init2(temp, mpwp);
  mpc_init2(g, mpwp);

  k = 0;
  for (i = 0; i < nclust; i++)
    k = MAX(k, punt[i + 1] - punt[i]);

  for (i = 0; i < nclust; i++) {	/* loop1: */
    if ((punt[i + 1] - punt[i]) == 1)
      continue;
    tst = true;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* looptst: */
      l = clust[punt[i] + j];
      if (!again[l])
	goto loop1;
      if (goal[0] == 'c') {
	if (status[l][0] == 'c' && status[l][2] == 'u') {
	  tst = false;
	  break;
	}
      } else if ((status[l][0] == 'c' && status[l][2] == 'u')
		 || (status[l][0] == 'c' && status[l][2] == 'i')) {
	tst = false;
	break;
      }
    }				/* for */

    if (tst)
      goto loop1;

    /* Compute super center sc and super radius sr */
    mpf_set_ui(srmp, 0);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      mpf_set_rdpe(rea, drad[l]);
      mpf_add(srmp, srmp, rea);
    }
    mpc_set_ui(sc, 0, 0);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      mpf_set_rdpe(rea, drad[l]);
      mpc_mul_f(temp, mroot[l], rea);
      mpc_add_eq(sc, temp);
    }
    mpc_div_eq_f(sc, srmp);
    rdpe_set(sr, rdpe_zero);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      mpc_sub(temp, sc, mroot[l]);
      mpc_get_cdpe(tmp, temp);
      cdpe_mod(rtmp, tmp);
      rdpe_add_eq(rtmp, drad[l]);
      if (rdpe_lt(sr, rtmp))
	rdpe_set(sr, rtmp);
    }
    
    if(DOLOG) {
      fprintf(logstr,"    MRESTART: clust=%d\n      sc=",i);
      mpc_out_str(logstr, 10, 10, sc);
      fprintf(logstr,"\n      sr=");
      rdpe_outln_str(logstr, sr);
    }
    
    /* Check the relative width of the cluster
     * If it is greater than 1 do not shift
     * and set status[:1)='c' that means
     * keep iterating Aberth's step. 
     * Check also the Newton-isolation of the cluster */

    mpc_get_cdpe(tmp, sc);
    cdpe_mod(rtmp, tmp);
    
    if(DOLOG){
      rdpe_div(rtmp2,sr,rtmp);
      fprintf(logstr,"      relative width=");
      rdpe_outln_str(logstr, rtmp2);
    }
    
    if (rdpe_gt(sr, rtmp)) {
      for (j = punt[i]; j < punt[i + 1]; j++)
	status[clust[j]][0] = 'c';
      if (DOLOG)
	fprintf(logstr, "    MRESTART: cluster %d relat. large: skip to the next component\n",i);
      goto loop1;
    }
    
    /* Now check the Newton isolation of the cluster */
    rdpe_set(rtmp2, rdpe_zero);
    for (k = 0; k < nclust; k++){
      if (k != i)
	for (j = 0; j < punt[k + 1] - punt[k]; j++) {
	  mpc_sub(temp, sc, mroot[clust[punt[k] + j]]);
	  mpc_get_cdpe(tmp, temp);
	  cdpe_mod(rtmp, tmp);
          rdpe_sub_eq(rtmp,drad[clust[punt[k] + j]]);
          rdpe_sub_eq(rtmp,sr);
	  rdpe_inv_eq(rtmp);
	  rdpe_add_eq(rtmp2,rtmp);
	}
    }
    rdpe_mul_eq(rtmp2,sr);
    rdpe_set_d(rtmp1, 0.3);

    if (rdpe_gt(rtmp2, rtmp1)) {
      for (jj = punt[i]; jj < punt[i + 1]; jj++)
	status[clust[jj]][0] = 'c';
      if (DOLOG) {
	fprintf(logstr, "    MRESTART: Cluster not Newton isolated:");
	fprintf(logstr, "              skip to the next component\n");
      }
      goto loop1;
    }
 
    if(DOLOG){
      fprintf(logstr,"    MRESTART: Approximations of cluster %d\n", i);
      for (j = 0; j < punt[i + 1] - punt[i]; j++) {
	l = clust[punt[i] + j];
	mpc_get_cdpe(tmp, mroot[l]);
	cdpe_out_str(logstr, tmp);
	fprintf(logstr,"  rad=");
        rdpe_outln_str(logstr,drad[l]);
      }
    }


    /* Compute the coefficients of the derivative of p(x) having order
     * equal to the multiplicity of the cluster -1. */

    for (j = 0; j <= n; j++)
      mpc_set(mfpc1[j], mfpc[j]);
    for (j = 1; j < punt[i + 1] - punt[i]; j++) {
      for (k = 0; k <= n - j; k++)
	mpc_mul_ui(mfpc1[k], mfpc1[k + 1], k + 1);
    }
    for (j = 0; j < n - (punt[i + 1] - punt[i]) + 2; j++) {
      mpc_get_cdpe(tmp, mfpc1[j]);
      cdpe_mod(dap1[j], tmp);
    }

    /* create the vectors needed if the polynomial is sparse */

    if (data_type[0] == 's') {
      for (j = 0; j < n - (punt[i + 1] - punt[i]) + 2; j++) {
	if (rdpe_ne(dap1[j], rdpe_zero))
	  spar1[j] = true;
	else
	  spar1[j] = false;
      }
      for (j = 0; j < n - (punt[i + 1] - punt[i]) + 1; j++)
	mpc_mul_ui(mfppc1[j], mfpc1[j + 1], j + 1);
    }
    /* Apply at most max_newt_it steps of Newton's iterations
     * to the above derivative starting from the super center
     * of the cluster. */
    mpc_set(g, sc);

    if (DOLOG) {
      fprintf(logstr, "    MRESTART: g before newton=");
      mpc_outln_str(logstr, 10, 30, g);
    }
    for (j = 0; j < max_newt_it; j++) {		/* loop_newt: */
      rdpe_set(rad, rdpe_zero);
      mnewton(n - (punt[i + 1] - punt[i]) + 1, g, rad, corr, mfpc1,
	      mfppc1, dap1, spar1, &cont);
      if (cont) {
	mpc_sub_eq(g, corr);
	if (DOLOG) {
	  fprintf(logstr, "    MRESTART: radius=");
	  rdpe_outln_str(logstr, rad);
	  fprintf(logstr, "    MRESTART: at iteration %d, g=", j);
	  mpc_outln_str(logstr, 10, 100, g);
	}
      } else
	break;
    }
    if (DOLOG)
      fprintf(logstr, "    MRESTART: performed %d Newton iter\n", j);
    if (j == max_newt_it) {
      if (DOLOG)
	fprintf(logstr, "Exceeded maximum Newton iterations in mrestart\n");
        goto loop1;
    }
    mpc_sub(temp, sc, g);
    mpc_get_cdpe(tmp, temp);
    cdpe_mod(rtmp, tmp);
    if (rdpe_gt(rtmp, sr)) {
      if (DOLOG)
	fprintf(logstr, "The gravity center falls outside the cluster\n");
      goto loop1;
    }
    
    /* shift the variable and compute new approximations */
    
    if (DOLOG)
      fprintf(logstr, "      MRESTART: call mshift\n");
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      mpc_get_cdpe(droot[l], mroot[l]);
    }
/*#D perform shift only if the new computed sr is smaller than old*0.25 */
    rdpe_mul_d(rtmp, sr, 0.25);
/*#D AGO99 Factors: 0.1 (MPS2.0), 0.5 (GIUGN98) */
    mshift(punt[i + 1] - punt[i], i, sr, g);
    if(rdpe_lt(sr, rtmp)){ /* Perform shift only if the new clust is smaller */
      mpc_get_cdpe(tmp, g);
      cdpe_mod(rtmp, tmp);
      rdpe_mul_eq(rtmp, mp_epsilon);
      rdpe_mul_eq_d(rtmp, 2);
      
      for (j = 0; j < punt[i + 1] - punt[i]; j++) {
	l = clust[punt[i] + j];
	mpc_set_cdpe(mroot[l], droot[l]);
	mpc_add_eq(mroot[l], g);
	cdpe_mod(rtmp1, droot[l]);
	rdpe_mul_d(drad[l], rtmp1, 2.0 * (punt[i + 1] - punt[i]));
	if (rdpe_lt(drad[l], rtmp))	
	  rdpe_set(drad[l], rtmp);
      }
    } else { 

     if(DOLOG) {
       fprintf(logstr,"    MRESTART: DO NOT PERFORM RESTART\n");
       fprintf(logstr,"    MRESTART: new radius of the cluster is larger\n");
     }
     
     goto loop1;
    }
  loop1:;
  }

  mpc_clear(g);
  tmpc_clear(temp);
  tmpc_clear(corr);
  tmpc_clear(sc);
  tmpf_clear(srmp);
  tmpf_clear(rea);
}

/**************************************************************
*                      SUBROUTINE FSHIFT                      *
***************************************************************
 This routine computes the first m+1 coefficients of the shifted 
 polynomial p(x+g), by performing m+1 Horner divisions.
 Then it computes the new starting approximations for the i-th
 cluster for i=i_clust by ing fstart and by updating root.
 The status vector is changed into 'o' for the components that
 belong to a cluster with relative radius less than eps.
 The status vector is changed into 'x' for the components that
 cannot be represented as float.
 **************************************************************/
void
fshift(int m, int i_clust, double clust_rad, cplx_t g, rdpe_t eps)
{
  int i, j;
  double prec, ag;
  cplx_t s;

  /* Perform divisions */

  prec = DBL_EPSILON;
  ag = cplx_mod(g);
  for (i = 0; i <= n; i++)
    cplx_set(fppc1[i], fpc[i]);
  for (i = 0; i <= m; i++) {
    cplx_set(s, fppc1[n]);
    for (j = n - 1; j >= i; j--) {
      cplx_mul_eq(s, g);
      cplx_add_eq(s, fppc1[j]);
      cplx_set(fppc1[j], s);
    }
    cplx_set(fppc[i], s);
  }

  /* start */
  for (i = 0; i <= m; i++)
    fap1[i] = cplx_mod(fppc[i]);

  fstart(m, i_clust, clust_rad, ag, eps, fap1);
}

/***********************************************************
*                   SUBROUTINE DSHIFT                      *
***********************************************************/
void
dshift(int m, int i_clust, rdpe_t clust_rad,
       cdpe_t g, rdpe_t eps)
{
  int i, j;
  rdpe_t prec, ag;
  cdpe_t s;

  rdpe_set_d(prec, DBL_EPSILON);
  cdpe_mod(ag, g);
  for (i = 0; i <= n; i++)
    cdpe_set(dpc1[i], dpc[i]);
  for (i = 0; i <= m; i++) {
    cdpe_set(s, dpc1[n]);
    for (j = n - 1; j >= i; j--) {
      cdpe_mul_eq(s, g);
      cdpe_add_eq(s, dpc1[j]);
      cdpe_set(dpc1[j], s);
    }
    cdpe_set(dpc2[i], s);
  }

  /* start */
  for (i = 0; i <= m; i++)
    cdpe_mod(dap1[i], dpc2[i]);

  dstart(m, i_clust, clust_rad, ag, eps, dap1);
}

/*******************************************************
*              SUBROUTINE MSHIFT                       *
*******************************************************/
void
mshift(int m, int i_clust, rdpe_t clust_rad, mpc_t g)
{
  int i, j, k;
  long int mpwp_temp, mpwp_max;
  rdpe_t ag, ap, abp, as, mp_ep;
  cdpe_t abd;
  mpc_t s;

  mpc_init2(s, mpwp);

  /* Perform divisions
   * In the mp version of the shift stage the computation
   * is performed with increasing levels of working precision
   * until the coefficients of the shifted polynomial have at
   * least one correct bit. */

  rdpe_set(mp_ep, mp_epsilon);
  mpc_get_cdpe(abd, g);
  cdpe_mod(ag, abd);
  for (i = 0; i <= n; i++)
    mpc_set(mfpc1[i], mfpc[i]);
  rdpe_set(as, rdpe_zero);
  rdpe_set(ap, rdpe_one);
  mpc_set_ui(s, 0, 0);
  k = 0;

  /* store the current working precision mpnw into mpnw_tmp */
  mpwp_temp = mpwp;
  mpwp_max = mpwp;

  do {				/* loop */
    mpc_set(s, mfpc1[n]);
    mpc_get_cdpe(abd, mfpc[n]);
    cdpe_mod(ap, abd);
    for (j = n - 1; j >= 0; j--) {
      mpc_get_cdpe(abd, mfpc[j]);
      cdpe_mod(abp, abd);
      rdpe_mul_eq(ap, ag);
      rdpe_mul_eq_d(abp, (double) j);
      rdpe_add_eq(ap, abp);
      mpc_mul_eq(s, g);
      mpc_add_eq(s, mfpc1[j]);
      mpc_set(mfpc1[j], s);
    }

    mpc_set(mfppc1[0], s);
    mpc_get_cdpe(abd, s);
    cdpe_mod(as, abd);
    rdpe_mul_eq(ap, mp_ep);
    rdpe_mul_eq_d(ap, 4.0 * (n + 1));
    k++;

    if (rdpe_lt(as, ap)) {
      mpwp_temp += mpwp;

      if (mpwp_temp > mpwp_max || mpwp_temp > prec_out * m * 2) { 
	if (DOLOG)
	  fprintf(logstr, "Reached the maximum allowed precision in mshift\n");
	break;
      }
      rdpe_set_2dl(mp_ep, 1.0, 1 - mpwp_temp);
      raisetemp(mpwp_temp);
      mpc_set_prec(s, (unsigned long int) mpwp_temp);
      mpc_set_prec(g, (unsigned long int) mpwp_temp);
      if (mpwp_max < mpwp_temp) 
	mpwp_max = mpwp_temp;

      for (j = 0; j <= n; j++)
	mpc_set(mfpc1[j], mfpc[j]);
    }
  } while (rdpe_lt(as, ap) && (k <= m));	/* loop */

  for (i = 1; i <= m; i++) {
    mpwp_temp = MAX(mpwp_temp - mpwp, mpwp);
    raisetemp_raw(mpwp_temp);
    mpc_set_prec_raw(s, (unsigned long int) mpwp_temp);
    mpc_set_prec_raw(g, (unsigned long int) mpwp_temp);
    mpc_set(s, mfpc1[n]);

    for (j = n - 1; j >= i; j--) {
      mpc_mul_eq(s, g);
      mpc_add_eq(s, mfpc1[j]);
      mpc_set(mfpc1[j], s);
    }
    mpc_set(mfppc1[i], s);

  }
  /*
    raisetemp_raw(mpwp);
    mpc_set_prec_raw(s, (unsigned long int) mpwp);
    mpc_set_prec_raw(g, (unsigned long int) mpwp);
  
   segue alternativa
  */
  raisetemp_raw(mpwp_max);
  mpc_set_prec_raw(s, (unsigned long int) mpwp_max);
  mpc_set_prec_raw(g, (unsigned long int) mpwp_max);
  raisetemp(mpwp);
  mpc_set_prec(s, (unsigned long int) mpwp);
  mpc_set_prec(g, (unsigned long int) mpwp);

  if (rdpe_lt(as, ap)) {
    for (j = 0; j < m; j++)
      rdpe_set(dap1[j], ap);
    mpc_get_cdpe(abd, mfppc1[m]);
    cdpe_mod(dap1[m], abd);
  } else
    for (i = 0; i <= m; i++) {
      mpc_get_cdpe(abd, mfppc1[i]);
      cdpe_mod(dap1[i], abd);
    }

  mstart(m, i_clust, clust_rad, ag, dap1);

  mpc_clear(s);
}

/**************************************************************
 *               SUBROUTINE RAISETEMP                         *
 *************************************************************/
void
raisetemp(unsigned long int digits)
{
  int i;

  for (i = 0; i <= n; i++) {
    mpc_set_prec(mfpc1[i], digits);
    mpc_set_prec(mfppc1[i], digits);
  }
}

/**************************************************************
 *               SUBROUTINE RAISETEMP_RAW                     *
 *************************************************************/
void
raisetemp_raw(unsigned long int digits)
{
  int i;

  for (i = 0; i <= n; i++) {
    mpc_set_prec_raw(mfpc1[i], digits);
    mpc_set_prec_raw(mfppc1[i], digits);
  }
}

/**************************************************************
 *               SUBROUTINE MNEWTIS                           *
 *************************************************************/
void 
mnewtis(void)
{
  boolean tst;
  int i, j, k, l, jj;
  rdpe_t sr, rtmp, rtmp1;
  cdpe_t tmp;
  tmpf_t rea, srmp;
  tmpc_t sc, temp;
  rdpe_t rtmp2; 

  /* For user's polynomials skip the restart stage (not yet implemented) */
  if (data_type[0] == 'u')
    return;
  tmpf_init2(rea, mpwp);
  tmpf_init2(srmp, mpwp);
  tmpc_init2(sc, mpwp);
  tmpc_init2(temp, mpwp);

  k = 0;
  for (i = 0; i < nclust; i++)
    k = MAX(k, punt[i + 1] - punt[i]);

  for (i = 0; i < nclust; i++) {	/* loop1: */

    if ((punt[i + 1] - punt[i]) == 1)
      continue;
    tst = true;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* looptst: */
      l = clust[punt[i] + j];
      if (!again[l])
	goto loop1;
      if (goal[0] == 'c') {
	if (status[l][0] == 'c' && status[l][2] == 'u') {
	  tst = false;
	  break;
	}
      } else if ((status[l][0] == 'c' && status[l][2] == 'u')
		 || (status[l][0] == 'c' && status[l][2] == 'i')) {
	tst = false;
	break;
      }
    }				/* for */
    if (tst)
      goto loop1;

    /* Compute super center sc and super radius sr */
    mpf_set_ui(srmp, 0);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      mpf_set_rdpe(rea, drad[l]);
      mpf_add(srmp, srmp, rea);
    }
    mpc_set_ui(sc, 0, 0);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      mpf_set_rdpe(rea, drad[l]);
      mpc_mul_f(temp, mroot[l], rea);
      mpc_add_eq(sc, temp);
    }
    mpc_div_eq_f(sc, srmp);
    rdpe_set(sr, rdpe_zero);
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      mpc_sub(temp, sc, mroot[l]);
      mpc_get_cdpe(tmp, temp);
      cdpe_mod(rtmp, tmp);
      rdpe_add_eq(rtmp, drad[l]);
      if (rdpe_lt(sr, rtmp))
	rdpe_set(sr, rtmp);
    }

    /* Check the relative width of the cluster
     * If it is greater than 1 do not shift
     * and set status[:1)='c' that means
     * keep iterating Aberth's step. 
     * Check also the Newton-isolation of the cluster */

    mpc_get_cdpe(tmp, sc);
    cdpe_mod(rtmp, tmp);
    rdpe_div(rtmp2,sr,rtmp);
    if (rdpe_gt(sr, rtmp)) {
      for (j = punt[i]; j < punt[i + 1]; j++)
	status[clust[j]][0] = 'c';
      if (DOLOG)
	fprintf(logstr, "   MNEWTIS cluster %d relat. large: "
		"skip to the next component\n", i);
      goto loop1;
    }
    
    /* Now check the Newton isolation of the cluster */
    rdpe_set(rtmp2, rdpe_zero);        
    for (k = 0; k < nclust; k++){
      if (k != i)
	for (j = 0; j < punt[k + 1] - punt[k]; j++) {
	  mpc_sub(temp, sc, mroot[clust[punt[k] + j]]);
	  mpc_get_cdpe(tmp, temp);
	  cdpe_mod(rtmp, tmp);
          rdpe_sub_eq(rtmp,drad[clust[punt[k] + j]]);
          rdpe_sub_eq(rtmp,sr);
	  rdpe_inv_eq(rtmp);
	  rdpe_add_eq(rtmp2,rtmp);
	}
    }
    rdpe_mul_eq(rtmp2,sr);
    rdpe_set_d(rtmp1,0.3);

    if (rdpe_gt(rtmp2, rtmp1)) {
      for (jj = punt[i]; jj < punt[i + 1]; jj++)
	status[clust[jj]][0] = 'c';
      if (DOLOG) {
	fprintf(logstr, "   MNEWTIS Cluster not Newton isolated:");
	fprintf(logstr, "           skip to the next component\n");
      }
      goto loop1;
    } 
    newtis=1;

  loop1:;
}
  
  tmpc_clear(temp);
  tmpc_clear(sc);
  tmpf_clear(srmp);
  tmpf_clear(rea);
}
