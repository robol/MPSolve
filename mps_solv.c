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
 * @brief Implementation of the routines to solve the polynomial
 * equation
 */

#include "mps.h"

/**
 @brief Set <code>again[i]</code> to <code>true</code> or to <code>false</code> 
 according to the values of <code>status[i,*]</code> and <code>goal</code>.

 More precisely:

 Goal count: .true. only for statu='**u' but not for 
         'f*u', 'a*u', 'o*u'                            (1) 
      - multipl. on: true also for 'c**'
      - Real on    : true also for '*u*' excluded (1)
      - Imag. on   : true also for '*v*' excluded (1)

 Goal isolate: true only for statu='**u' or statu='c**' but
          not for 'f**', 'a**', 'o**', 'i*i', 'i*o'     (2) 
      - multipl. on: true also for 'c**'
      - Real on    : true also for '*u*' excluded (2)
      - Imag. on   : true also for '*v*' excluded (2)
 Goal approximate: the same as isolate.
*/
void
update(void)
{
  int i;

  for (i = 0; i < n; i++)
    again[i] = false;
  switch (goal[0]) {

  case 'c':			/*  count */
    for (i = 0; i < n; i++) {
      if (status[i][2] == 'u')
	if (status[i][0] != 'f' && status[i][0] != 'a'
	    && status[i][0] != 'o')
	  again[i] = true;
      if (goal[2] == 'm' && status[i][0] == 'c' && status[i][0] != 'o')
	again[i] = true;

      switch (goal[3]) {

      case 'r':		/* real option */
	if (status[i][1] == 'w'
	    && (status[i][2] != 'u'
		|| (status[i][0] != 'f' && status[i][0] != 'a'
		    && status[i][0] != 'o')))
	  again[i] = true;
	break;

      case 'i':		/* imaginary option */
	if (status[i][1] == 'w'
	    && (status[i][2] != 'u'
		|| (status[i][0] != 'f' && status[i][0] != 'a'
		    && status[i][0] != 'o')))
	  again[i] = true;
	break;

      case 'b':		/* both imaginary and real options */
	if (status[i][1] == 'w'
	    && (status[i][2] != 'u'
		|| (status[i][0] != 'f' && status[i][0] != 'a'
		    && status[i][0] != 'o')))
	  again[i] = true;
	break;
      }
    }

    break;

  case 'i':			/* isolate */
    for (i = 0; i < n; i++) {
      if (status[i][2] == 'u' || (status[i][0] == 'c' && status[i][2] == 'i'))
	if (status[i][0] != 'f' && status[i][0] != 'a' && status[i][0] != 'o'
	    && (status[i][0] != 'i' || status[i][2] != 'i'))
	  again[i] = true;
      if (goal[2] == 'm' && status[i][0] == 'c' && status[i][2] != 'o')
	again[i] = true;

      switch (goal[3]) {

      case 'r':		/* real option */
	if (status[i][1] == 'w'
	    && (status[i][0] != 'f' && status[i][0] != 'a'
		&& status[i][0] != 'o'))
	  again[i] = true;
	break;

      case 'i':		/* imaginary option */
	if (status[i][1] == 'w' &&
	    (status[i][0] != 'f' && status[i][0] != 'a' &&
	     status[i][0] != 'o'))
	  again[i] = true;	/* DARIO RIVEDERE */
	break;

      case 'b':		/* both imaginary and real options */
	if (status[i][1] == 'w'
	    && (status[i][0] != 'f'
		&& status[i][0] != 'a' && status[i][0] != 'o'))
	  again[i] = true;
	break;
      }
    }

    break;

  case 'a':			/* approximate (the same as isolate) */
    for (i = 0; i < n; i++) {
      if (status[i][2] == 'u' || (status[i][0] == 'c' && status[i][2] == 'i'))
	if (status[i][0] != 'f' && status[i][0] != 'a' && status[i][0] != 'o')
	  again[i] = true;

      if (goal[2] == 'm' && status[i][0] == 'c' && status[i][2] != 'o')
	again[i] = true;

      switch (goal[3]) {

      case 'r':		/* real option */
	if (status[i][1] == 'w'
	    && (status[i][0] != 'f' && status[i][0] != 'a'
		&& status[i][0] != 'o'))
	  again[i] = true;
	break;

      case 'i':		/* imaginary option */
	if (status[i][1] == 'w'
	    && (status[i][0] != 'f'
		&& status[i][0] != 'a' && status[i][0] != 'o'))
	  again[i] = true;
	break;

      case 'b':		/* both imaginary and real options */
	if (status[i][1] == 'w'
	    && (status[i][0] != 'f'
		&& status[i][0] != 'a' && status[i][0] != 'o'))
	  again[i] = true;
	break;
      }
    }
    break;
  }
}

/********************************************************
*                    SUBROUTINE FSRAD                   *
*********************************************************
 Compute super center and super radius: float case
 Compute the super radius of the i-th cluster
*********************************************************/
void
fsrad(int i, cplx_t sc, double *sr)
{
  cplx_t ctmp;
  double sum;
  int j, l;

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
  *sr = 0.0;
  for (j = 0; j < punt[i + 1] - punt[i]; j++) {
    l = clust[punt[i] + j];
    cplx_sub(ctmp, sc, froot[l]);
    *sr = MAX(*sr, frad[l] + cplx_mod(ctmp));
  }
}

/********************************************************
*                    SUBROUTINE DSRAD                   *
*********************************************************
 Compute super center and super radius: DPE case
 Compute the super radius of the i-th cluster
*********************************************************/
void
dsrad(int i, cdpe_t sc, rdpe_t sr)
{
  cdpe_t ctmp;
  rdpe_t sum, rtmp;
  int j, l;

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
}

/********************************************************
*                    SUBROUTINE MSRAD                   *
*********************************************************
 Compute super center and super radius: MP case
 Compute the super radius of the i-th cluster
*********************************************************/
void
msrad(int i, mpc_t sc, rdpe_t sr)
{
  int j, l;
  rdpe_t rtmp;
  cdpe_t cdtmp;
  tmpf_t ftmp, sum;
  tmpc_t ctmp;

  tmpc_init2(ctmp, mpwp);
  tmpf_init2(ftmp, mpwp);
  tmpf_init2(sum, mpwp);

  mpf_set_ui(sum, 0);
  for (j = 0; j < punt[i + 1] - punt[i]; j++) {
    l = clust[punt[i] + j];
    mpf_set_rdpe(ftmp, drad[l]);
    mpf_add(sum, sum, ftmp);
  }

  mpc_set_ui(sc, 0, 0);
  for (j = 0; j < punt[i + 1] - punt[i]; j++) {
    l = clust[punt[i] + j];
    mpf_set_rdpe(ftmp, drad[l]);
    mpc_mul_f(ctmp, mroot[l], ftmp);
    mpc_add_eq(sc, ctmp);
  }

  mpc_div_eq_f(sc, sum);
  rdpe_set(sr, rdpe_zero);
  for (j = 0; j < punt[i + 1] - punt[i]; j++) {
    l = clust[punt[i] + j];
    mpc_sub(ctmp, sc, mroot[l]);
    mpc_get_cdpe(cdtmp, ctmp);
    cdpe_mod(rtmp, cdtmp);
    rdpe_add_eq(rtmp, drad[l]);
    if (rdpe_lt(sr, rtmp))
      rdpe_set(sr, rtmp);
  }

  tmpf_clear(sum);
  tmpf_clear(ftmp);
  tmpc_clear(ctmp);
}

/******************************************************
*           SUBROUTINE FMODIFY                        *
*******************************************************
 Modify the vector 'statu' according to the goal, and
 to the location of the roots.
 
 The subroutine is used also for marking the new cluster
 that have been detected between two consecutive packets
 of Aberth's iteration.
 ==1== The subroutine changes into 'C' the components of
 status[:1) corresponding to old clusters, keeping
 status[:,1)='c' for the new formed clusters.
 In this way applying restart selects new starting
 approximations only for the new detected clusters.
 ==2== For the components for which 
 status[:,1]!='C', 'f', 'x' performs the following
 analysis:
 If the cluster has mult=1 mark it with status[:1)='i'
    if is also approximated mark it with status(:1)='a'
 Check if c*u and i*u (i.e., uncertain set) can 
 be made certain according to goal[1]
 
 ==3== 
 Perform the same with options, that is,
 If multiplicity is on then check if a cluster
   corresponds to a  multiple root
 If detect real then detect real roots
 if detect imaginary then detect imaginary roots
 If detect both then detect both imaginary and
   real roots
**************************************************/
void
fmodify(void)
{
  int i, j, l, k, nnewclust, i_new, i_old, s, ip1, i1, l1, j1, nf, j2,
   l2;
  double sr, tmpr, afri, sep1;
  cplx_t sc;
  boolean tcr, tcr1;

  /* ==1== Change into 'C' the components of status for old clusters */
  nf = 2 * n;  /* Isolation factor */
  for (i = 0; i < n; i++)
    if (status[i][0] == 'c')
      status[i][0] = 'C';

  i_old = 0;
  i_new = 0;
  s = 1;
  for (i = 1; i <= nclust && i_new < n; i++) {	/*  loop1: DO i=1, nclust */
    if (oldpunt[i_old + 1] == punt[i_new + 1]) {
      i_old++;
      i_new++;
      continue;
    } else {
      for (j = i_new + 1; j < nclust; j++) {	/* loop2: DO j=i_new+1, nclust */
	if (oldpunt[i_old + 1] != punt[j + 1])
	  continue;
	else {
	  nnewclust = j - i_new + 1;	/* scan each new cluster */
	  for (k = 0; k < nnewclust; k++) {	/* loop3: DO k=1, nnewclust */
	    i1 = i_new + k;
      /********************************* 
         scan the entries of each new cluster 
         set status[l][0]='i' if the cluster has multip=1
         and mark with 'c'
         the ones which are different from 'i'
       **********************************/
	    if (punt[i1 + 1] - punt[i1] == 1 && status[clust[punt[i1]]][0] != 'x'
		&& status[clust[punt[i1]]][0] != 'f')
	      status[clust[punt[i1]]][0] = 'i';
	    for (l = 0; l < punt[i1 + 1] - punt[i1]; l++) {	/* loop4: */
	      ip1 = clust[punt[i1] + l];
	      if (status[ip1][0] != 'i' && status[ip1][0] != 'x'
		  && status[ip1][0] != 'f' &&
		  status[ip1][0] != 'a' &&
		  status[ip1][0] != 'o')
		status[ip1][0] = 'c';
	    }
	  }
	  i_new = j + 1;
	  i_old++;
	  break;
	}			/* else */
      }				/* for */
    }				/* else */
  }				/* for */

/*=2== Scan all the clusters */
  for (i = 0; i < nclust; i++) {	/* scan : DO i=1,nclust */
    /* check isolation/approximation */
    if (punt[i + 1] - punt[i] == 1 && status[clust[punt[i]]][0] != 'x'
	&& status[clust[punt[i]]][0] != 'f') {
      status[clust[punt[i]]][0] = 'i';
      tmpr = cplx_mod(froot[clust[punt[i]]]);
      tmpr = frad[clust[punt[i]]] / tmpr;
      tmpr = log(tmpr);
      if (tmpr < -prec_out * LOG2)
	status[clust[punt[i]]][0] = 'a';
    }
    /* Scan inside the cluster */
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan_in */
      l = clust[punt[i] + j];
      /**************************************************
  first check for inside/outside unit disk in the case where
  there are very large and/or very small roots (statu='x')
  and for counting only
  *************************************************/
      afri = cplx_mod(froot[i]);
      if (status[l][0] == 'x' && goal[0] == 'c') {
	if ((goal[1] == 'i' && afri < 1) ||
	    (goal[1] == 'o' && afri > 1))
	  status[l][2] = 'i';
	if ((goal[1] == 'i' && afri > 1) ||
	    (goal[1] == 'o' && afri < 1))
	  status[l][2] = 'o';
      }
      /* Now check the standard cases */
      if (status[l][0] == 'x')
	continue;
      if ((status[l][0] == 'c' || status[l][0] == 'i'
	   || status[l][0] == 'C' || status[l][0] == 'o')
	  &&status[l][2] == 'u') {
	/* Check if the approximation is inside/outside the set */
	switch (goal[1]) {

	case 'a':		/* all */
	  status[l][2] = 'i';
	  break;

	case 'i':		/* inside unit circle */
	  if (!ftouchunit(nf, l)) {
	    if (cplx_mod(froot[l]) < 1)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'o':		/* outside unit circle */
	  if (!ftouchunit(nf, l)) {
	    if (cplx_mod(froot[l]) > 1)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'l':		/* left half plane  */
	  if (!ftouchimag(nf, l)) {
	    if (cplx_Re(froot[l]) < 0)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'r':		/* right half plane */
	  if (!ftouchimag(nf, l)) {
	    if (cplx_Re(froot[l]) > 0)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'u':		/* upper half plane */
	  if (!ftouchreal(nf, l)) {
	    if (cplx_Im(froot[l]) > 0)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'd':		/* lower half plane */
	  if (!ftouchreal(nf, l)) {
	    if (cplx_Im(froot[l]) < 0)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'R':		/* Real line  NEW */
	  if (status[l][1] != 'w')
	    continue;
	  if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	    if (ftouchreal(1, l)) {
	      if (data_type[1] == 'r') {
		status[l][2] = 'i';
		status[l][1] = 'R';
	      } else {
		/* fsrad(i, sc, &sr); DARIO */
		sr = frad[l];
		if (log(sr) < sep - n * lmax_coeff) {
		  status[l][2] = 'i';
		  status[l][1] = 'R';
		} else {
		  status[l][2] = 'u';
		  status[l][1] = 'w';
		}
	      }
	    } else {		/* do not touch real */
	      status[l][2] = 'o';
	      status[l][1] = 'r';
	    }
	    continue;
	  } else {		/* cluster */
	    tcr = ftouchreal(nf, l);
	    for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	      l1 = clust[punt[i] + j1];
	      tcr1 = ftouchreal(nf, l1);
	      if ((tcr && tcr1) || (!tcr && !tcr1))
		continue;
	      else {		/*  mixed situation */
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
		goto scan;
	      }
	    }
	    if (tcr) {		/* tutti i dischi intersecano R */
	      if (data_type[2] == 'f' || data_type[2] == 'b') {
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
	      } else {		/* integer/rational polynomial */
		fsrad(i, sc, &sr);
		sep1 = sep;
		if (data_type[1] == 'c')
		  sep1 = sep - n * lmax_coeff;
		if (log(sr) < sep1) {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'i';
		    status[l2][1] = 'R';
		  }
		} else {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'u';
		    status[l2][1] = 'w';
		  }
		}
	      }
	    } else {		/* tutti i dischi non intersecano R */

		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'o';
		  status[l2][1] = 'r';
		}
	    }
	  }
	  break;

	case 'I':		/* Imaginary line  NEW */
	  if (status[l][1] != 'w' && status[l][1] != 'r')
	    continue;
	  if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	    if (ftouchimag(nf, l)) {
	      /* fsrad(i, sc, &sr); DARIO */
	      sr = frad[l];
	      if (log(sr) < sep - n * lmax_coeff) {
		status[l][2] = 'i';
		status[l][1] = 'I';
	      } else {
		status[l][2] = 'u';
		status[l][1] = 'w';
	      }
	    } else {		/* do not touch imag */
	      status[l][2] = 'o';
	      status[l][1] = 'i';
	    }
	    continue;
	  } else {		/* cluster */
	    tcr = ftouchimag(nf, l);
	    for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	      l1 = clust[punt[i] + j1];
	      tcr1 = ftouchimag(nf, l1);
	      if ((tcr && tcr1) || (!tcr && !tcr1))
		continue;
	      else {		/* mixed situation */
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
		goto scan;
	      }
	    }
	    if (tcr) {		/* tutti i dischi intersecano I */
	      if (data_type[2] == 'f' || data_type[2] == 'b') {
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
	      } else {		/* integer/rational polynomial */
		fsrad(i, sc, &sr);
		sep1 = sep;
		sep1 = sep - n * lmax_coeff;
		if (log(sr) < sep1) {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'i';
		    status[l2][1] = 'I';
		  }
		} else {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'u';
		    status[l2][1] = 'w';
		  }
		}
	      }
	    } else {		/* tutti i dischi non intersecano I */
	      status[l][2] = 'o';
	      status[l][1] = 'i';
	    }
	  }
	  break;

	case 'S':		/* Set provided by the user */
	  error(1, "User set not implemented yet");
	  break;

	default:
	  error(1, "Mistake in goal");
	  break;
	}
      }
    }
    /* If some cluster still contains an uncertain disk then set
     * all the disks uncertain */
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      if (status[l][2] == 'u') {
	for (j1 = 0; j1 < punt[i + 1] - punt[i]; j1++) {
	  l1 = clust[punt[i] + j1];
	  status[l1][2] = 'u';
	}
	break;
      }
    }
  scan:;			/* scan the next component */
  }

/*==3==  now check the options */

  /* Option multiplicity */
  for (i = 0; i < nclust; i++) {	/* scan1 */
    if (punt[i + 1] - punt[i] == 1)
      continue;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan1_in */
      l = clust[punt[i] + j];
      if (status[l][0] == 'x' ||
	  status[l][0] == 'f' || status[l][0] == 'm')
	goto scan1;
      if (goal[2] == 'm' && (status[l][0] == 'c' ||	/* NEW */
			     status[l][0] == 'C')) {	/* multiplicity on */
	if (data_type[2] == 'b' || data_type[2] == 'f')	/* float coeff. */
	  error(1, "Fatal: Float coefficients - impossible to detect multiplicity");

	/* compute super center and super radius */
	fsrad(i, sc, &sr);

	if (log(sr) < sep)
	  for (j1 = 0; j1 < punt[i + 1] - punt[i]; j1++) {
	    l1 = clust[punt[i] + j1];	/* NEW j-> j1 */
	    status[l1][0] = 'm';
	  }
	goto scan1;
      }
    }
  scan1:;			/* scan next component */
  }

  /* Option Real check */
  if ((goal[3] == 'r' || goal[3] == 'b') && goal[1] != 'R') {
    for (i = 0; i < nclust; i++) {	/* scan2 */
      for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan2_in */
	l = clust[punt[i] + j];
	if (status[l][1] != 'w' && status[l][1] != 'i')
	  continue;
	if (status[l][0] == 'x' || status[l][0] == 'f')
	  goto scan2;
	if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	  if (ftouchreal(nf, l)) {
	    if (data_type[1] == 'r')
	      status[l][1] = 'R';
	    else {
	      /* fsrad(i, sc, &sr); DARIO */
	      sr = frad[l];
	      if (log(sr) < sep - n * lmax_coeff)
		status[l][1] = 'R';
	      else
		status[l][1] = 'w';
	    }
	  } else		/* do not touch real */
	    status[l][1] = 'r';
	  continue;
	} else {		/* cluster */
	  tcr = ftouchreal(nf, l);
	  for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	    l1 = clust[punt[i] + j1];
	    tcr1 = ftouchreal(nf, l1);
	    if ((tcr && tcr1) || (!tcr && !tcr1))
	      continue;
	    else {		/* mixed situation */
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	      }
	      goto scan2_in;
	    }
	  }
	  if (tcr) {		/* tutti i dischi intersecano R */
	    if (data_type[2] == 'f' || data_type[2] == 'b')
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	    } else {		/* integer/rational polynomial */
	      fsrad(i, sc, &sr);
	      sep1 = sep;
	      if (data_type[1] == 'c')
		sep1 = sep - n * lmax_coeff;
	      if (log(sr) < sep1)
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'R';
	      } else
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'w';
		}
	    }
	  } else		/* tutti i dischi non intersecano R */
	    status[l][1] = 'r';
	}
      scan2_in:;
      }
    scan2:;
    }
  }
  /* Option Imaginary check */
  if (goal[3] == 'i' || goal[3] == 'b') {
    for (i = 0; i < nclust; i++) {	/* scan3 */
      for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan3_in */
	l = clust[punt[i] + j];
	if (status[l][0] == 'x' || status[l][0] == 'f')
	  goto scan3;
	if (status[l][1] != 'w' && status[l][1] != 'r')
	  continue;
	if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	  if (ftouchimag(nf, l)) {
	    /* fsrad(i, sc, &sr); DARIO */
	    sr = frad[l];
	    if (log(sr) < sep - n * lmax_coeff)
	      status[l][1] = 'I';
	    else
	      status[l][1] = 'w';
	  } else		/* do not touch imag */
	    status[l][1] = 'i';
	  continue;
	} else {		/* cluster */
	  tcr = ftouchimag(nf, l);
	  for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	    l1 = clust[punt[i] + j1];
	    tcr1 = ftouchimag(nf, l1);
	    if ((tcr && tcr1) || (!tcr && !tcr1))
	      continue;
	    else {		/* mixed situation */
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	      }
	      goto scan3_in;
	    }
	  }
	  if (tcr) {		/* tutti i dischi intersecano I */
	    if (data_type[2] == 'f' || data_type[2] == 'b')
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	    } else {		/* integer/rational polynomial */
	      fsrad(i, sc, &sr);
	      sep1 = sep;
	      sep1 = sep - n * lmax_coeff;
	      if (log(sr) < sep1)
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'I';
	      } else
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'w';
		}
	    }
	  } else		/* tutti i dischi non intersecano I */
	    status[l][1] = 'i';
	}
      scan3_in:;
      }
    scan3:;
    }
  }
}

/****************************************************
*           SUBROUTINE DMODIFY                      *
****************************************************/
void
dmodify(void)
{
  int i, j, l, k, nnewclust, i_new, i_old, s, ip1, i1, l1, j1, j2, l2,
   nf;
  double rtmp, sep1;
  rdpe_t sr, tmpr;
  cdpe_t sc;
  boolean tcr, tcr1;

/* ==1==  Change into 'C' the components of status for old clusters */
  nf = 2 * n; /* Isolation factor */
  for (i = 0; i < n; i++)
    if (status[i][0] == 'c')
      status[i][0] = 'C';

  i_old = 0;
  i_new = 0;
  s = 1;
  for (i = 1; i <= nclust && i_new < n; i++) {	/*  loop1: DO i=1, nclust */
    if (oldpunt[i_old + 1] == punt[i_new + 1]) {
      i_old++;
      i_new++;
      continue;
    } else {
      for (j = i_new + 1; j < nclust; j++) {	/* loop2:  */
	if (oldpunt[i_old + 1] != punt[j + 1])
	  continue;
	else {
	  nnewclust = j - i_new + 1;	/*  scan each new cluster */
	  for (k = 0; k < nnewclust; k++) {	/* loop3:  */
	    i1 = i_new + k;
	    /* scan the entries of each new cluster 
	     * set status[l][0]='i' if the cluster has multip=1
	     * and mark with 'c' those which are different from 'i'
	     */
	    if (punt[i1 + 1] - punt[i1] == 1 && status[clust[punt[i1]]][0] != 'x'
		&& status[clust[punt[i1]]][0] != 'f')
	      status[clust[punt[i1]]][0] = 'i';
	    for (l = 0; l < punt[i1 + 1] - punt[i1]; l++) {	/* loop4: */
	      ip1 = clust[punt[i1] + l];
	      if (status[ip1][0] != 'i' && status[ip1][0] != 'x'
		  && status[ip1][0] != 'f' &&
		  status[ip1][0] != 'a' &&
		  status[ip1][0] != 'o')
		status[ip1][0] = 'c';
	    }
	  }
	  i_new = j + 1;
	  i_old++;
	  break;
	}
      }
    }
  }

/*=2== Scan all the clusters */
  for (i = 0; i < nclust; i++) {	/* scan: */
    /* check isolation/approximation */
    if (punt[i + 1] - punt[i] == 1 && status[clust[punt[i]]][0] != 'x'
	&& status[clust[punt[i]]][0] != 'f') {
      status[clust[punt[i]]][0] = 'i';
      cdpe_mod(tmpr, droot[clust[punt[i]]]);
      rdpe_div(tmpr, drad[clust[punt[i]]], tmpr);
      rtmp = rdpe_log(tmpr);
      if (rtmp < -prec_out * LOG2)
	status[clust[punt[i]]][0] = 'a';
    }
    /* Scan inside the cluster */
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan_in: */
      l = clust[punt[i] + j];

      /* Now check the standard cases */
      if (status[l][0] == 'x')
	continue;
      if ((status[l][0] == 'c' || status[l][0] == 'i'
	   || status[l][0] == 'C' || status[l][0] == 'o')
	  &&status[l][2] == 'u') {
	/* Check if the approximation is inside/outside the set */
	switch (goal[1]) {

	case 'a':		/* all */
	  status[l][2] = 'i';
	  break;

	case 'i':		/* inside unit circle */
	  if (!dtouchunit(nf, l)) {
	    cdpe_mod(tmpr, droot[l]);
	    if (rdpe_lt(tmpr, rdpe_one))
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'o':		/* outside unit circle */
	  if (!dtouchunit(nf, l)) {
	    cdpe_mod(tmpr, droot[l]);
	    if (rdpe_gt(tmpr, rdpe_one))
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'l':		/* left half plane  */
	  if (!dtouchimag(nf, l)) {
	    rdpe_set(tmpr, cdpe_Re(droot[l]));
	    if (rdpe_lt(tmpr, rdpe_zero))
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'r':		/* right half plane */
	  if (!dtouchimag(nf, l)) {
	    rdpe_set(tmpr, cdpe_Re(droot[l]));
	    if (rdpe_gt(tmpr, rdpe_zero))
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'u':		/* upper half plane */
	  if (!dtouchreal(nf, l)) {
	    rdpe_set(tmpr, cdpe_Im(droot[l]));
	    if (rdpe_gt(tmpr, rdpe_zero))
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'd':		/* lower half plane */
	  if (!dtouchreal(nf, l)) {
	    rdpe_set(tmpr, cdpe_Im(droot[l]));
	    if (rdpe_lt(tmpr, rdpe_zero))
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'R':		/* Real line  NEW */
	  if (status[l][1] != 'w')
	    continue;
	  if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	    if (dtouchreal(1, l)) {
	      if (data_type[1] == 'r') {
		status[l][2] = 'i';
		status[l][1] = 'R';
	      } else {
		/* dsrad(i, sc, sr); DARIO */
		rdpe_set(sr, drad[l]);
		if (rdpe_log(sr) < sep - n * lmax_coeff) {
		  status[l][2] = 'i';
		  status[l][1] = 'R';
		} else {
		  status[l][2] = 'u';
		  status[l][1] = 'w';
		}
	      }
	    } else {		/* do not touch real */
	      status[l][2] = 'o';
	      status[l][1] = 'r';
	    }
	    continue;
	  } else {		/* cluster */
	    tcr = dtouchreal(nf, l);
	    for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	      l1 = clust[punt[i] + j1];
	      tcr1 = dtouchreal(nf, l1);
	      if ((tcr && tcr1) || (!tcr && !tcr1))
		continue;
	      else {		/*  mixed situation */
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
		goto scan;
	      }
	    }
	    if (tcr) {		/* tutti i dischi intersecano R */
	      if (data_type[2] == 'f' || data_type[2] == 'b') {
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
	      } else {		/* integer/rational polynomial */
		dsrad(i, sc, sr);
		sep1 = sep;
		if (data_type[1] == 'c')
		  sep1 = sep - n * lmax_coeff;
		if (rdpe_log(sr) < sep1) {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'i';
		    status[l2][1] = 'R';
		  }
		} else {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'u';
		    status[l2][1] = 'w';
		  }
		}
	      }
	    } else {		/* tutti i dischi non intersecano R */

		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'o';
		  status[l2][1] = 'r';
		}
	    }
	  }
	  break;

	case 'I':		/* Imaginary line  NEW */
	  if (status[l][1] != 'w' && status[l][1] != 'r')
	    continue;
	  if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	    if (dtouchimag(nf, l)) {
	      /* dsrad(i, sc, sr); */
	      rdpe_set(sr, drad[l]);
	      if (rdpe_log(sr) < sep - n * lmax_coeff) {
		status[l][2] = 'i';
		status[l][1] = 'I';
	      } else {
		status[l][2] = 'u';
		status[l][1] = 'w';
	      }
	    } else {		/* do not touch imag */
	      status[l][2] = 'o';
	      status[l][1] = 'i';
	    }
	    continue;
	  } else {		/* cluster */
	    tcr = dtouchimag(nf, l);
	    for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	      l1 = clust[punt[i] + j1];
	      tcr1 = dtouchimag(nf, l1);
	      if ((tcr && tcr1) || (!tcr && !tcr1))
		continue;
	      else {		/* mixed situation */
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
		goto scan;
	      }
	    }
	    if (tcr) {		/* tutti i dischi intersecano I */
	      if (data_type[2] == 'f' || data_type[2] == 'b') {
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
	      } else {		/* integer/rational polynomial */
		dsrad(i, sc, sr);
		sep1 = sep;
		sep1 = sep - n * lmax_coeff;
		if (rdpe_log(sr) < sep1) {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'i';
		    status[l2][1] = 'I';
		  }
		} else {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'u';
		    status[l2][1] = 'w';
		  }
		}
	      }
	    } else {		/* tutti i dischi non intersecano I */
	      status[l][2] = 'o';
	      status[l][1] = 'i';
	    }
	  }
	  break;

	case 'S':		/* Set provided by the user */
	  error(1, "Custom region not implemented yet");
	  break;
	default:
	  error(1, "mistake in goal");
	  break;
	}
      }
    }
    /* If some cluster still contains an uncertain disk then set
     * all the disks uncertain
     */
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      if (status[l][2] == 'u') {
	for (j1 = 0; j1 < punt[i + 1] - punt[i]; j1++) {
	  l1 = clust[punt[i] + j1];
	  status[l1][2] = 'u';
	}
	break;
      }
    }
  scan:;
  }

/*==3==  now check the options */

  /* Option multiplicity */
  for (i = 0; i < nclust; i++) {	/* scan1 */
    if (punt[i + 1] - punt[i] == 1)
      continue;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan1_in */
      l = clust[punt[i] + j];
      if (status[l][0] == 'x' ||
	  status[l][0] == 'f' || status[l][0] == 'm')
	goto scan1;
      if (goal[2] == 'm' && (status[l][0] == 'c' ||	/* NEW */
			     status[l][0] == 'C')) {	/* multiplicity on */

	if (data_type[2] == 'b' || data_type[2] == 'f') /* float coeff. */
	  error(1, "Fatal: Float coefficients - impossible to detect multiplicity");

	/* compute super center and super radius */
	dsrad(i, sc, sr);

	if (rdpe_log(sr) < sep)
	  for (j1 = 0; j1 < punt[i + 1] - punt[i]; j1++) {
	    l1 = clust[punt[i] + j1];	/* NEW j-> j1 */
	    status[l1][0] = 'm';
	  }
	goto scan1;
      }
    }
  scan1:;			/* scan next component */
  }

  /* Option Real check */
  if ((goal[3] == 'r' || goal[3] == 'b') && goal[1] != 'R') {
    for (i = 0; i < nclust; i++) {	/* scan2 */
      for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan2_in */
	l = clust[punt[i] + j];
	if (status[l][1] != 'w' && status[l][1] != 'i')
	  continue;
	if (status[l][0] == 'x' || status[l][0] == 'f')
	  goto scan2;
	if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	  if (dtouchreal(nf, l)) {
	    if (data_type[1] == 'r')
	      status[l][1] = 'R';
	    else {
	      /* dsrad(i, sc, sr); DARIO */
	      rdpe_set(sr, drad[l]);
	      if (rdpe_log(sr) < sep - n * lmax_coeff)
		status[l][1] = 'R';
	      else
		status[l][1] = 'w';
	    }
	  } else		/* do not touch real */
	    status[l][1] = 'r';
	  continue;
	} else {		/* cluster */
	  tcr = dtouchreal(nf, l);
	  for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	    l1 = clust[punt[i] + j1];
	    tcr1 = dtouchreal(nf, l1);
	    if ((tcr && tcr1) || (!tcr && !tcr1))
	      continue;
	    else {		/* mixed situation */
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	      }
	      goto scan2_in;
	    }
	  }
	  if (tcr) {		/* tutti i dischi intersecano R */
	    if (data_type[2] == 'f' || data_type[2] == 'b')
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	    } else {		/* integer/rational polynomial */
	      dsrad(i, sc, sr);
	      sep1 = sep;
	      if (data_type[1] == 'c')
		sep1 = sep - n * lmax_coeff;
	      if (rdpe_log(sr) < sep1)
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'R';
	      } else
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'w';
		}
	    }
	  } else		/* tutti i dischi non intersecano R */
	    status[l][1] = 'r';
	}
      scan2_in:;
      }
    scan2:;
    }
  }

  /* Option Imaginary check */
  if (goal[3] == 'i' || goal[3] == 'b') {
    for (i = 0; i < nclust; i++) {	/* scan3 */
      for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan3_in */
	l = clust[punt[i] + j];
	if (status[l][0] == 'x' || status[l][0] == 'f')
	  goto scan3;
	if (status[l][1] != 'w' && status[l][1] != 'r')
	  continue;
	if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	  if (dtouchimag(nf, l)) {
	    /* dsrad(i, sc, sr); DARIO */
	    rdpe_set(sr, drad[l]);
	    if (rdpe_log(sr) < sep - n * lmax_coeff)
	      status[l][1] = 'I';
	    else
	      status[l][1] = 'w';
	  } else		/* do not touch imag */
	    status[l][1] = 'i';
	  continue;
	} else {		/* cluster */
	  tcr = dtouchimag(nf, l);
	  for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	    l1 = clust[punt[i] + j1];
	    tcr1 = dtouchimag(nf, l1);
	    if ((tcr && tcr1) || (!tcr && !tcr1))
	      continue;
	    else {		/* mixed situation */
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	      }
	      goto scan3_in;
	    }
	  }
	  if (tcr) {		/* tutti i dischi intersecano I */
	    if (data_type[2] == 'f' || data_type[2] == 'b')
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	    } else {		/* integer/rational polynomial */
	      dsrad(i, sc, sr);
	      sep1 = sep;
	      sep1 = sep - n * lmax_coeff;
	      if (rdpe_log(sr) < sep1)
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'I';
	      } else
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'w';
		}
	    }
	  } else		/* tutti i dischi non intersecano I */
	    status[l][1] = 'i';
	}
      scan3_in:;
      }
    scan3:;
    }
  }
}

/****************************************************
*           SUBROUTINE MMODIFY                      *
****************************************************/
void
mmodify(void)
{
  int i, j, l, k, nnewclust, i_new, i_old, s, ip1, i1, l1, j1, nf, j2,
   l2;
  double rtmp, sep1;
  rdpe_t sr, tmpr;
  cdpe_t tmpc;
  boolean tcr, tcr1;
  tmpf_t tmpf;
  tmpc_t sc;

  tmpc_init2(sc, mpwp);
  tmpf_init2(tmpf, mpwp);

  /* ==1==  Change into 'C' the components of status for old clusters */
  nf = 2 * n; /* Isolation factor */
  for (i = 0; i < n; i++)
    if (status[i][0] == 'c')
      status[i][0] = 'C';

  i_old = 0;
  i_new = 0;
  s = 1;
  for (i = 1; i <= nclust && i_new < n; i++) {	/* loop1: */
    if (oldpunt[i_old + 1] == punt[i_new + 1]) {
      i_old++;
      i_new++;
      continue;
    } else {
      for (j = i_new + 1; j < nclust; j++) {	/* loop2: */
	if (oldpunt[i_old + 1] != punt[j + 1])
	  continue;
	else {
	  nnewclust = j - i_new + 1;	/* scan each new cluster */
	  for (k = 0; k < nnewclust; k++) {	/* loop3: */
	    i1 = i_new + k;
      /*****************************************
        scan the entries of each new cluster set
        status[l][0]='i' if the cluster has multip=1 and
        mark with 'c' the ones which are different from 'i' 
      *****************************************/
	    if (punt[i1 + 1] - punt[i1] == 1 && status[clust[punt[i1]]][0] != 'x'
		&& status[clust[punt[i1]]][0] != 'f')
	      status[clust[punt[i1]]][0] = 'i';
	    for (l = 0; l < punt[i1 + 1] - punt[i1]; l++) {	/* loop4: */
	      ip1 = clust[punt[i1] + l];
	      if (status[ip1][0] != 'i' && status[ip1][0] != 'x'
		  && status[ip1][0] != 'f' &&
		  status[ip1][0] != 'a' &&
		  status[ip1][0] != 'o')
		status[ip1][0] = 'c';
	    }
	  }
	  i_new = j + 1;
	  i_old++;
	  break;
	}
      }
    }
  }

/*=2== Scan all the clusters */
  for (i = 0; i < nclust; i++) {	/*  scan : DO i=1,nclust */
    /* check isolation/approximation */
    if (punt[i + 1] - punt[i] == 1 && status[clust[punt[i]]][0] != 'x'
	&& status[clust[punt[i]]][0] != 'f') {
      status[clust[punt[i]]][0] = 'i';
      mpc_get_cdpe(tmpc, mroot[clust[punt[i]]]);
      cdpe_mod(tmpr, tmpc);
      rdpe_div(tmpr, drad[clust[punt[i]]], tmpr);
      rtmp = rdpe_log(tmpr);
      if (rtmp < -prec_out * LOG2)
	status[clust[punt[i]]][0] = 'a';
    }
    /* Scan inside the cluster */
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan_in: */
      l = clust[punt[i] + j];

      /* Now check the standard cases */
      if (status[l][0] == 'x')
	continue;
      if ((status[l][0] == 'c' || status[l][0] == 'i' || status[l][0] == 'a'
	   || status[l][0] == 'C' || status[l][0] == 'o')
	  &&status[l][2] == 'u') {
	/* Check if the approximation is inside/outside the set */
	switch (goal[1]) {

	case 'a':		/* all */
	  status[l][2] = 'i';
	  break;

	case 'i':		/* inside unit circle */
	  if (!mtouchunit(nf, l)) {
	    mpc_get_cdpe(tmpc, mroot[l]);
	    cdpe_mod(tmpr, tmpc);
	    if (rdpe_lt(tmpr, rdpe_one))
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'o':		/* outside unit circle */
	  if (!mtouchunit(nf, l)) {
	    mpc_get_cdpe(tmpc, mroot[l]);
	    cdpe_mod(tmpr, tmpc);
	    if (rdpe_gt(tmpr, rdpe_one))
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'l':		/* left half plane  */
	  if (!mtouchimag(nf, l)) {
	    /* mpc_get_re(tmpf, mroot[l]); */
	    mpf_set(tmpf, mpc_Re(mroot[l]));
	    if (mpf_sgn(tmpf) == -1)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'r':		/* right half plane */
	  if (!mtouchimag(nf, l)) {
	    /* mpc_get_re(tmpf, mroot[l]); */
	    mpf_set(tmpf, mpc_Re(mroot[l]));
	    if (mpf_sgn(tmpf) == 1)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'u':		/* upper half plane */
	  if (!mtouchreal(nf, l)) {
	    /* mpc_get_im(tmpf, mroot[l]); */
	    mpf_set(tmpf, mpc_Im(mroot[l]));
	    if (mpf_sgn(tmpf) == 1)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'd':		/* lower half plane */
	  if (!mtouchreal(nf, l)) {
	    /* mpc_get_im(tmpf, mroot[l]); */
	    mpf_set(tmpf, mpc_Im(mroot[l]));
	    if (mpf_sgn(tmpf) == -1)
	      status[l][2] = 'i';
	    else
	      status[l][2] = 'o';
	  }
	  break;

	case 'R':		/* Real line  NEW */
	  if (status[l][1] != 'w')
	    continue;
	  if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	    if (mtouchreal(1, l)) {
	      if (data_type[1] == 'r') {
		status[l][2] = 'i';
		status[l][1] = 'R';
	      } else {
		rdpe_set(sr, drad[l]);
		/* msrad(i, sc, sr);*/ /*#DARIO*/
		if (rdpe_log(sr) < sep - n * lmax_coeff) {
		  status[l][2] = 'i';
		  status[l][1] = 'R';
		} else {
		  status[l][2] = 'u';
		  status[l][1] = 'w';
		}
	      }
	    } else {		/* do not touch real */
	      status[l][2] = 'o';
	      status[l][1] = 'r';
	    }
	    continue;
	  } else {		/* cluster */
	    tcr = mtouchreal(nf, l);
	    for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	      l1 = clust[punt[i] + j1];
	      tcr1 = mtouchreal(nf, l1);
	      if ((tcr && tcr1) || (!tcr && !tcr1))
		continue;
	      else {		/*  mixed situation */
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
		goto scan;
	      }
	    }
	    if (tcr) {		/* tutti i dischi intersecano R */
	      if (data_type[2] == 'f' || data_type[2] == 'b') {
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
	      } else {		/* integer/rational polynomial */
		msrad(i, sc, sr);
		sep1 = sep;
		if (data_type[1] == 'c')
		  sep1 = sep - n * lmax_coeff;
		if (rdpe_log(sr) < sep1) {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'i';
		    status[l2][1] = 'R';
		  }
		} else {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'u';
		    status[l2][1] = 'w';
		  }
		}
	      }
	    } else {		/* tutti i dischi non intersecano R */

		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'o';
		  status[l2][1] = 'r';
		}
	    }
	  }
	  break;

	case 'I':
	  if (status[l][1] != 'w' && status[l][1] != 'r')
	    continue;
	  if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	    if (mtouchimag(nf, l)) {
	      /* msrad(i, sc, sr);*/ /*#DARIO */
	      rdpe_set(sr, drad[l]);
	      if (rdpe_log(sr) < sep - n * lmax_coeff) {
		status[l][2] = 'i';
		status[l][1] = 'I';
	      } else {
		status[l][2] = 'u';
		status[l][1] = 'w';
	      }
	    } else {		/* do not touch imag */
	      status[l][2] = 'o';
	      status[l][1] = 'i';
	    }
	    continue;
	  } else {		/* cluster */
	    tcr = mtouchimag(nf, l);
	    for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	      l1 = clust[punt[i] + j1];
	      tcr1 = mtouchimag(nf, l1);
	      if ((tcr && tcr1) || (!tcr && !tcr1))
		continue;
	      else {		/* mixed situation */
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
		goto scan;
	      }
	    }
	    if (tcr) {		/* tutti i dischi intersecano I */
	      if (data_type[2] == 'f' || data_type[2] == 'b') {
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][2] = 'u';
		  status[l2][1] = 'w';
		}
	      } else {		/* integer/rational polynomial */
		msrad(i, sc, sr);
		sep1 = sep;
		sep1 = sep - n * lmax_coeff;
		if (rdpe_log(sr) < sep1) {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'i';
		    status[l2][1] = 'I';
		  }
		} else {
		  for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		    l2 = clust[punt[i] + j2];
		    status[l2][2] = 'u';
		    status[l2][1] = 'w';
		  }
		}
	      }
	    } else {		/* tutti i dischi non intersecano I */
	      status[l][2] = 'o';
	      status[l][1] = 'i';
	    }
	  }
	  break;

	case 'S':		/* Set provided by the user */
	  error(1, "Custom region not implemented yet");
	  break;
	default:
	  error(1, "mistake in goal");
	  break;
	}
      }
    }
    /* If some cluster still contains an uncertain disk then set
     * all the disks uncertain
     */
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {
      l = clust[punt[i] + j];
      if (status[l][2] == 'u') {
	for (j1 = 0; j1 < punt[i + 1] - punt[i]; j1++) {
	  l1 = clust[punt[i] + j1];
	  status[l1][2] = 'u';
	}
	break;
      }
    }

  scan:;
  }

/*==3==  now check the options */

  /* Option multiplicity */
  for (i = 0; i < nclust; i++) {	/* scan1 */
    if (punt[i + 1] - punt[i] == 1)
      continue;
    for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan1_in */
      l = clust[punt[i] + j];
      if (status[l][0] == 'x' ||
	  status[l][0] == 'f' || status[l][0] == 'm')
	goto scan1;
      if (goal[2] == 'm' && (status[l][0] == 'c' ||	/* NEW */
			     status[l][0] == 'C')) {	/* multiplicity on */
			     
	if (data_type[2] == 'b' || data_type[2] == 'f') /* float coeff. */
	  error(1, "Fatal: Float coefficients - impossible to detect multiplicity");

	/* compute super center and super radius */
	msrad(i, sc, sr);

	if (rdpe_log(sr) < sep)
	  for (j1 = 0; j1 < punt[i + 1] - punt[i]; j1++) {
	    l1 = clust[punt[i] + j1];	/* NEW j-> j1 */
	    status[l1][0] = 'm';
	  }
	goto scan1;
      }
    }
  scan1:;			/* scan next component */
  }

  /* Option Real check */
  if ((goal[3] == 'r' || goal[3] == 'b') && goal[1] != 'R') {
    for (i = 0; i < nclust; i++) {	/* scan2 */
      for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan2_in */
	l = clust[punt[i] + j];
	if (status[l][1] != 'w' && status[l][1] != 'i')
	  continue;
	if (status[l][0] == 'x' || status[l][0] == 'f')
	  goto scan2;
	if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	  if (mtouchreal(nf, l)) {
	    if (data_type[1] == 'r')
	      status[l][1] = 'R';
	    else {
	      /* msrad(i, sc, sr);*/ /*#DARIO */
	      rdpe_set(sr, drad[l]);
	      if (rdpe_log(sr) < sep - n * lmax_coeff)
		status[l][1] = 'R';
	      else
		status[l][1] = 'w';
	    }
	  } else		/* do not touch real */
	    status[l][1] = 'r';
	  continue;
	} else {		/* cluster */
	  tcr = mtouchreal(nf, l);
	  for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	    l1 = clust[punt[i] + j1];
	    tcr1 = mtouchreal(nf, l1);
	    if ((tcr && tcr1) || (!tcr && !tcr1))
	      continue;
	    else {		/* mixed situation */
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	      }
	      goto scan2_in;
	    }
	  }
	  if (tcr) {		/* tutti i dischi intersecano R */
	    if (data_type[2] == 'f' || data_type[2] == 'b')
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	    } else {		/* integer/rational polynomial */
	      msrad(i, sc, sr);
	      sep1 = sep;
	      if (data_type[1] == 'c')
		sep1 = sep - n * lmax_coeff;
	      if (rdpe_log(sr) < sep1)
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'R';
	      } else
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'w';
		}
	    }
	  } else		/* tutti i dischi non intersecano R */
	    status[l][1] = 'r';
	}
      scan2_in:;
      }
    scan2:;
    }
  }
  /* Option Imaginary check */
  if (goal[3] == 'i' || goal[3] == 'b') {
    for (i = 0; i < nclust; i++) {	/* scan3 */
      for (j = 0; j < punt[i + 1] - punt[i]; j++) {	/* scan3_in */
	l = clust[punt[i] + j];
	if (status[l][0] == 'x' || status[l][0] == 'f')
	  goto scan3;
	if (status[l][1] != 'w' && status[l][1] != 'r')
	  continue;
	if (punt[i + 1] - punt[i] == 1) {	/* one disk */
	  if (mtouchimag(nf, l)) {
	    rdpe_set(sr, drad[l]);
	    /* msrad(i, sc, sr); DARIO */
	    if (rdpe_log(sr) < sep - n * lmax_coeff)
	      status[l][1] = 'I';
	    else
	      status[l][1] = 'w';
	  } else		/* do not touch imag */
	    status[l][1] = 'i';
	  continue;
	} else {		/* cluster */
	  tcr = mtouchimag(nf, l);
	  for (j1 = 1; j1 < punt[i + 1] - punt[i]; j1++) {
	    l1 = clust[punt[i] + j1];
	    tcr1 = mtouchimag(nf, l1);
	    if ((tcr && tcr1) || (!tcr && !tcr1))
	      continue;
	    else {		/* mixed situation */
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	      }
	      goto scan3_in;
	    }
	  }
	  if (tcr) {		/* tutti i dischi intersecano I */
	    if (data_type[2] == 'f' || data_type[2] == 'b')
	      for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		l2 = clust[punt[i] + j2];
		status[l2][1] = 'w';
	    } else {		/* integer/rational polynomial */
	      msrad(i, sc, sr);
	      sep1 = sep;
	      sep1 = sep - n * lmax_coeff;
	      if (rdpe_log(sr) < sep1)
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'I';
	      } else
		for (j2 = 0; j2 < punt[i + 1] - punt[i]; j2++) {
		  l2 = clust[punt[i] + j2];
		  status[l2][1] = 'w';
		}
	    }
	  } else		/* tutti i dischi non intersecano I */
	    status[l][1] = 'i';
	}
      scan3_in:;
      }
    scan3:;
    }
  }

  tmpf_clear(tmpf);
  tmpc_clear(sc);
}

/*********************************************************
*           SUBROUTINE CHECK_STOP                        *
**********************************************************
 Set computed= true if the stop condition is satisfied,
 Set computed= false otherwise.
 The stop condition is obtained from the vector status as
 follows:
 If the Goal is count:
    == Stop if:
      --  **u does not exist, except for a*u, o*u, f*u
      -- Mult. and does not exist c**
      -- Real. and does not exist *u*, except for au*, ou*
      -- Imag  and does not exist *v*, except for av*, ov*
 If the Goal is isolate or approximate:
    == Stop if:
      -- **u does not exist, except for  a*u, o*u, f*u
         and if c*i does not exist.
      -- Mult. and does not exist c*i, o*i
      -- Real. and does not exist *ui, except for aui, oui
      -- Imag  and does not exist *vi, except for avi, ovi
************************************************************/
boolean
check_stop(void)
{
  int i;
  boolean computed;

  computed = false;
  /* count */
  if (goal[0] == 'c') {
    for (i = 0; i < n; i++) {
      if (status[i][2] == 'u' && status[i][0] != 'f' && status[i][0] != 'o'
	  && status[i][0] != 'a')
	return computed;
      if (goal[2] == 'm' && status[i][0] == 'c' && status[i][2] != 'o')
	return computed;
      if (goal[3] == 'r' && status[i][1] == 'w' && status[i][2] != 'o'
	  && status[i][0] != 'a' && status[i][0] != 'o'
	  && status[i][0] != 'm')	/* NEW */
	return computed;
      if (goal[3] == 'i' && status[i][1] == 'w' && status[i][2] != 'o'
       && status[i][0] != 'a' && status[i][0] != 'o' && status[i][0] != 'm')
	return computed;
      if (goal[3] == 'b' && status[i][2] != 'o'
	&& status[i][1] == 'w' && status[i][0] != 'a' && status[i][0] != 'o'
	  && status[i][0] != 'm')
	return computed;
    }
    computed = true;
  }
  /* isolate or approximate */
  if (goal[0] == 'i' || goal[0] == 'a') {
    for (i = 0; i < n; i++) {
      if (status[i][2] == 'u' && status[i][0] != 'f' && status[i][0] != 'o'
	  && status[i][0] != 'a')
	return computed;
      if (status[i][0] == 'c' && status[i][2] != 'o')
	return computed;
      if (goal[2] == 'm' && status[i][2] != 'o'
	  && (status[i][0] == 'c'))
	return computed;
      if (goal[3] == 'r' && status[i][1] == 'w' && status[i][2] != 'o'
       && status[i][0] != 'a' && status[i][0] != 'o' && status[i][0] != 'm')
	return computed;
      if (goal[3] == 'i' && status[i][1] == 'w' && status[i][2] != 'o' && status[i][0] != 'm'
       && status[i][0] != 'a' && status[i][0] != 'o' && status[i][0] != 'm')
	return computed;
      if (goal[3] == 'b' && status[i][2] != 'o' && status[i][0] != 'm'
       && status[i][1] == 'w' && status[i][0] != 'a' && status[i][0] != 'o')
	return computed;
    }
    computed = true;
  }
  return computed;
}

/**
  @brief Actually solve the polynomial

 This routine performs the following computations:
 -# Select starting approximations and check if some of them need dpe.
    Initialize the vector again which is true if the corresponding
    approximation is out of the root neighbourhood.
 -# Performs max_pack packets of Aberth iterations on all the
    components out of the root neighbourhood belonging to the set S
    and on which it is possible to iterate with float.
    More precisely, each packet performs max_it iterations on all the
    components where again is true. At each iteration check if the
    current approximation is in the root neighbourhood; in this case
    set 'again' to false.  
 -# If at the end of the general packet all the approximations
    are inside the root neighbourhood, i.e., 'again' is false in all
    the components then return.
    else, perform cluster analysis, select new starting approximations
    update the vector 'statu', update the vector 'again' that selects 
    the components on which to iterate, according to the goal, and 
    repeat until the max number of allowed packets is reached. 
    In the latter case output FAILURE.

 The local variable 'again' controls the iteration: i.e., 
   again[i]= true means iterate on the i-th component
*/
void 
fsolve(boolean * d_after_f)
{
  boolean excep;
  int it_pack, iter, nit, oldnclust, i, j;
  rdpe_t eps_out;

  /* == 1 ==  Initialize variables */
  it_pack = 0;
  nclust = 0;
  for (i = 0; i < n; i++) {
    again[i] = true;
    cplx_set(froot[i], cplx_zero);
    frad[i] = DBL_MAX;
  }

  /* choose starting approximations */
  if (DOLOG)
    fprintf(logstr, "FSOLVE: call fstart");

  fstart(n, 0, 0.0, 0.0, eps_out, fap);
  /***************
     this part of code performs shift in the gravity center of the roots 
     In order to use it, uncomment the part below and comment the 
     instruction above. Dangerous for overflow.
  ************/
  /*
  {
    cplx_t ft;
    cplx_mul_d(ft, fpc[n], -n);
    cplx_div(ft, fpc[n-1], ft);
    fshift(n, 0, 100, ft, eps_out);
  }
  */ /* till here */ 
  
  if (DOLOG)
    dump(logstr);

  /* Check if there are too large or too small approximations */
  *d_after_f = false;
  for (i = 0; i < n; i++)
    if (status[i][0] == 'x') {
      again[i] = false;
      *d_after_f = true;
    }

  /* == 2 ==  Perform max_pack packets of Aberth's iterations */
  if (DOLOG)
    fprintf(logstr, "   FSOLVE:  call fpolzer\n");
  for (iter = 0; iter < max_pack; iter++) {	/* floop: */

    fpolzer(&nit, &excep);
    it_pack += nit;

    if (DOLOG)
      fprintf(logstr, "Packet %d  iterations= %d\n", iter, nit);

    /* perform cluster analysis, shift, restart, update 'statu', and
     * update 'again'
     */
    if (excep) {
      for (i = 0; i <= n; i++)
	oldpunt[i] = punt[i];
      oldnclust = nclust;

      if (DOLOG)
	fprintf(logstr, "   FSOLVE: call fcluster\n");
      /* cluster analysis */
      fcluster(2 * n); /* Isolation factor */
      if (oldnclust == nclust) {
	if (DOLOG)
	  fprintf(logstr, "   FSOLVE: cycle\n");
	continue;
      } else {
	/* modify the vector status and mark also the old
	 * clusters with 'C'
	 */
	if (DOLOG)
	  fprintf(logstr, "   FSOLVE: call modify\n");
	fmodify();

	if (iter == 0)
	  for (i = 0; i < n; i++)
	    if (status[i][0] == 'C')
	      status[i][0] = 'c';

	/* If the polynomial is not given in terms of its coeff. then
	 * skip the restart stage */
	if (data_type[0] != 'u') {
	  /* choose new starting approximations only for new clusters */
	  if (DOLOG)
	    fprintf(logstr, "   FSOLVE: call frestart\n");
	  frestart();
	}
	/* reset the status vector */
	for (j = 0; j < n; j++) {
	  if (status[j][0] == 'C')
	    status[j][0] = 'c';
	  again_old[j] = again[j];
	}

	/* update 'again' */
	if (DOLOG)
	  fprintf(logstr, "   FSOLVE: call update\n");
	update();

	/* adjust 'again' This is needed since we are
	 * between two packets */
	for (i = 0; i < n; i++)
	  if (!again_old[i])
	    again[i] = false;
	if (DOLOG)
	  fprintf(logstr, "   FSOLVE: call checkstop\n");
	/* Check the stop condition */
	if (check_stop())
	  return;
      }
    } else
      break;
  }

  /* The 'floop' has been completed: 
   * If the max number of iteration has been reached then output FAILURE */
  if (iter == max_pack) {
    dump(logstr);
    error(1, "Float: reached the maximum number of packet iterations");
  }
  /* Otherwise exit since all the approximations are
   * in the root neighbourhood, except for the ones that cannot be
   * represented as double. */

  if (DOLOG)
    fprintf(logstr, "FLOAT: nit= %d\n", it_pack);

  /* Update */
  if (DOLOG)
    fprintf(logstr, "   FSOLVE: call fcluster\n");
  for (i = 0; i <= n; i++)
    oldpunt[i] = punt[i];
  oldnclust = nclust;
  fcluster(2 * n); /* Isolation factor */

  if (DOLOG)
    fprintf(logstr, "   FSOLVE: call modify\n");
  fmodify();

  /* reset the status vector */
  for (j = 0; j < n; j++)
    if (status[j][0] == 'C')
      status[j][0] = 'c';
}

/**********************************************************
*          SUBROUTINE FPOLZER                             *
***********************************************************
 This routine applies nit iterations of Aberth's method to
 the i-th component of the approximations for which again[i]
 is true. Set again[i]=false if the i-th approximation is in
 the root neighbourhood. Stop if  again[i]=false for any i.
 excep= true if after nit iterations some approximation is
 still out of the root neighbourhood.
************************************************************/
void
fpolzer(int *it, boolean * excep)
{
  int i, iter, nzeros;
  cplx_t corr, abcorr;
  double rad1, modcorr;

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  nzeros = 0;
  for (i = 0; i < n; i++)
    if (!again[i])
      nzeros++;
  if (nzeros == n)
    return;

  /* Start Aberth's iterations */
  if (DOLOG)
    fprintf(logstr, "FPOLZER: starts aberth it\n");

  for (iter = 0; iter < max_it; iter++) {	/* do_iter : DO iter=1,nit */

    if (DOLOG) {
      fprintf(logstr, "FPOLZER: iteration %d\n", iter);
      dump(logstr);
    }

    for (i = 0; i < n; i++) {	/* do_index */

      if (again[i]) {
	(*it)++;
	rad1 = frad[i];
	if (data_type[0] != 'u') {
	  fnewton(n, froot[i], &frad[i], corr, fpc, fap, &again[i]);
	  if (iter == 0 && !again[i] && frad[i] > rad1 && rad1 != 0)
	    frad[i] = rad1;
	  /***************************************
	  The above condition is needed to cope with the case
	  where at the first iteration the starting point
	  is already in the root neighbourhood and the actually
	  computed radius is too big since the value of the first
	  derivative is too small.
	  In this case the previous radius bound, obtained by
	  means of Rouche' is more reliable and strict
	  **************************************/
	} else
	  fnewton_usr(froot[i], &frad[i], corr, &again[i]);
	  
	if (again[i] ||
	  /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
	    data_type[0] == 'u' || iter != 0 || frad[i] != rad1) {
	  faberth(i, abcorr);
	  cplx_mul_eq(abcorr, corr);
	  cplx_sub(abcorr, cplx_one, abcorr);
	  cplx_div(abcorr, corr, abcorr);
          cplx_sub_eq(froot[i], abcorr);
	  modcorr = cplx_mod(abcorr);
	  frad[i] += modcorr;
	}

	/* check for new approximated roots */
	if (!again[i]) {
	  nzeros++;
	  if (nzeros == n)
	    return;
	}

      }
    }
  }
  *excep = true;
}

/**************************************************************
*              SUBROUTINE DPOLZER                             *
***************************************************************
 This routine applies nit iterations of Aberth's method to the
 i-th component of the approximations for which again[i] is true
 Set again[i]=false if the i-th approximation is in the root 
 neighbourhood
 Stop if  again[i]=false for any i.
 excep= true if after nit iterations some approximation is still
 out of the root neighbourhood.
***************************************************************/
void
dpolzer(int *it, boolean * excep)
{
  int iter, i, nzeros;
  rdpe_t rad1, rtmp;
  cdpe_t corr, abcorr;

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  nzeros = 0;
  for (i = 0; i < n; i++)
    if (!again[i])
      nzeros++;
  if (nzeros == n)
    return;

  /* Start Aberth's iterations */
  if (DOLOG)
    fprintf(logstr, "DPOLZER: starts aberth\n");
  for (iter = 0; iter < max_it; iter++) {	/* do_iter: */

    for (i = 0; i < n; i++) {	/* do_index: */

      if (again[i]) {
	(*it)++;
	rdpe_set(rad1, drad[i]);
	if (data_type[0] != 'u') {
	  dnewton(n, droot[i], drad[i], corr, dpc, dap, &again[i]);
	  if (iter == 0 && !again[i] && rdpe_gt(drad[i], rad1)
	      && rdpe_ne(rad1, rdpe_zero))
	    rdpe_set(drad[i], rad1);
	} else
	  dnewton_usr(droot[i], drad[i], corr, &again[i]);

  /************************************************
    The above condition is needed to manage with the case where
    at the first iteration the starting point is already in the
    root neighbourhood and the actually computed radius is too
    big since the value of the first derivative is too small.
    In this case the previous radius bound, obtained by means of
    Rouche' is more reliable and strict
    **********************************************/

	if (again[i] ||
	    /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
	    data_type[0] == 'u' || iter != 0 || rdpe_ne(drad[i], rad1)) {
	  daberth(i, abcorr);
	  cdpe_mul_eq(abcorr, corr);
	  cdpe_sub(abcorr, cdpe_one, abcorr);
	  cdpe_div(abcorr, corr, abcorr);
	  cdpe_sub_eq(droot[i], abcorr);
	  cdpe_mod(rtmp, abcorr);
	  rdpe_add_eq(drad[i], rtmp);
	}

	/* check for new approximated roots */
	if (!again[i]) {
	  nzeros++;
	  if (nzeros == n)
	    return;
	}

      }
    }
  }
  *excep = true;
}

/*****************************************************
*             SUBROUTINE DSOLVE                      *
******************************************************
 This routine applies nit iterations of Aberth's method
 to the i-th component of the approximations for which
 again[i] is true
 Set again[i]=false if the i-th approximation is in the
 root neighbourhood
 Stop if  again[i]=false for any i.
 excep= true if after nit iterations some approximation
 is still out of the root neighbourhood.
 ******************************************************/
void
dsolve(boolean d_after_f)
{
  int it_pack, iter, nit, oldnclust, i, j;
  boolean excep;
  rdpe_t dummy;

  if (DOLOG)
    fprintf(logstr, "   DSOLVE: d_after_f= %d\n", d_after_f);

  /* == 1 == Initialize variables */
  it_pack = 0;

  if (d_after_f)
    for (i = 0; i < n; i++)
      if (status[i][0] == 'x') {
	again[i] = true;
	rdpe_set_d(drad[i], DBL_MAX);
      } else
	again[i] = false;
  else {
    nclust = 0;
    for (i = 0; i < n; i++) {
      again[i] = true;
      rdpe_set(drad[i], RDPE_MAX);
      cdpe_set(droot[i], cdpe_zero);
    }
  }

  /* Choose starting approximations */
  if (DOLOG)
    fprintf(logstr, "   DSOLVE: call dstart con again=\n");

  rdpe_set(dummy, rdpe_zero);
  dstart(n, 0, dummy, dummy, dummy, dap);

  /* Now adjust the status vector */
  if (d_after_f)
    for (i = 0; i < n; i++)
      if (status[i][0] == 'x')
	status[i][0] = 'c';

  /* == 2 == Perform max_pack  packets of Aberth's iterations */
  if (DOLOG)
    fprintf(logstr, "   DSOLVE: call dpolzero\n");

  for (iter = 0; iter < max_pack; iter++) {	/* dloop : DO iter=1,max_pack */
    dpolzer(&nit, &excep);
    it_pack += nit;

    if (DOLOG)
      fprintf(logstr, "Packet %d iterations= %d\n", iter, nit);

    if (excep) {
      for (i = 0; i <= n; i++)
	oldpunt[i] = punt[i];
      oldnclust = nclust;

      /* cluster analysis */
      if (DOLOG)
	fprintf(logstr, "   DSOLVE: call dcluster\n");
      dcluster(2 * n); /* Isolation factor */
      if (oldnclust == nclust) {
	if (DOLOG)
	  fprintf(logstr, "   DSOLVE:  CYCLE\n");
	continue;
      } else {
	if (DOLOG)
	  fprintf(logstr, "   DSOLVE: call dmodify\n");
	dmodify();

	if (iter == 0 && !d_after_f)
	  for (i = 0; i < n; i++)
	    if (status[i][0] == 'C')
	      status[i][0] = 'c';

	/* If the polynomial is not given in terms of its
	 * coeff. then skip the restart stage */
	if (data_type[0] != 'u') {
	  /* choose new starting approximations only for new clusters */
	  if (DOLOG)
	    fprintf(logstr, "   DSOLVE: call drestart\n");
	  drestart();
	}
	/* reset the status vector */
	for (j = 0; j < n; j++)
	  if (status[j][0] == 'C')
	    status[j][0] = 'c';
	for (j = 0; j < n; j++)
	  again_old[j] = again[j];

	/* update 'again' */
	if (DOLOG)
	  fprintf(logstr, "   DSOLVE: call update\n");
	update();
	/* adjust 'again'
	 * This is needed since we are between two packets
	 */
	for (i = 0; i < n; i++)
	  if (!again_old[i])
	    again[i] = false;
	if (DOLOG)
	  fprintf(logstr, "   DSOLVE: call checkstop\n");
	if (check_stop())
	  return;
      }
    } else
      break;
  }
  if (iter == max_pack) {
    dump(logstr);
    error(1, "DPE: reached the maximum number of packet iterations");
  }
  /* Otherwise exit since all the approximations are
   * in the root neighbourhood, except for the ones that cannot be
   * represented as double. 
   */

  if (DOLOG)
    fprintf(logstr, "DPE: nit=%d\n", nit);

  /* Update */
  if (DOLOG)
    fprintf(logstr, "   DSOLVE: now update: call dcluster\n");
  for (i = 0; i <= n; i++)
    oldpunt[i] = punt[i];
  oldnclust = nclust;
  dcluster(2 * n); /* Isolation factor */

  if (DOLOG)
    fprintf(logstr, "   DSOLVE: now call dmodify\n");
  dmodify();

  /* reset the status vector */
  for (j = 0; j < n; j++)
    if (status[j][0] == 'C')
      status[j][0] = 'c';
}

/****************************************************
*             SUBROUTINE MSOLVE                     *
****************************************************/
void
msolve(void)
{
  int iter, nit, oldnclust, i, j, it_pack; 
  boolean excep;
  int nzc;

  /* == 1 == Initialize variables */
  it_pack = 0;

  if (DOLOG)
    fprintf(logstr, "  MSOLVE: call restart\n");
  if (data_type[0] != 'u')
    mrestart();
  if (DOLOG)
    fprintf(logstr, "  MSOLVE: call update1\n");
  update();
  if (DOLOG) {
    fprintf(logstr, "  MSOLVE: again after update = ");
    for (i = 0; i < n; i++)
      fprintf(logstr, "%d", again[i]);
    fprintf(logstr, "\n");
  }
  
  for (i = 0; i < n; i++)
    if (again[i])
      rootwp[i] = mpwp;

  if (DOLOG)
    fprintf(logstr, "  MSOLVE: call checkstop\n");
  if (check_stop()) {
    mmodify();

    /* reset the status vector */
    for (j = 0; j < n; j++)
      if (status[j][0] == 'C')
	status[j][0] = 'c';

    return;
  }
  nzc = 0;
  if (goal[1] == 'a' && goal[2] == 'n' && goal[3] == 'n')
    for (i = 0; i < n; i++)
      if (status[i][0] == 'i' || status[i][0] == 'a' || status[i][0] == 'o')
	nzc++;
  if (DOLOG)
    fprintf(logstr, "  MSOLVE: nzc=%d\n", nzc);

  if (nzc == n) {
    if (DOLOG)
      fprintf(logstr, "  MSOLVE: call mmodify and return\n");
    mmodify();

    /* reset the status vector */
    for (j = 0; j < n; j++)
      if (status[j][0] == 'C')
	status[j][0] = 'c';
    return;
  }

  /* Perform max_pack  packets of Aberth's iterations */
  if (DOLOG)
    fprintf(logstr, "  MSOLVE: Perform packets of Aberth\n");

  for (iter = 0; iter < max_pack; iter++) {	/* mloop : DO iter=1,max_pack */
    if (DOLOG) {
      fprintf(logstr, "  MSOLVE: packet= %d\n", iter);
      fprintf(logstr, "  MSOLVE: again before mpolzer =");
      for (i = 0; i < n; i++)
	fprintf(logstr, "%d", again[i]);
      fprintf(logstr, "\n");
      fprintf(logstr, "  MSOLVE: call mpolzer\n");
    }
    mpolzer(&nit, &excep);

    if (DOLOG)
      fprintf(logstr, "  MSOLVE: Packet %d: iterations= %d\n", iter, nit);

    it_pack += nit;
    nzc = 0;
    if (goal[1] == 'a' && goal[2] == 'n' && goal[3] == 'n')	/* DARIO APRILE 98 */
      for (i = 0; i < n; i++)
	if (status[i][0] == 'i' || status[i][0] == 'a' || status[i][0] == 'o')
	  nzc++;
    if (DOLOG)
      fprintf(logstr, "  MSOLVE: check again nzc=%d\n", nzc);
    if (nzc == n) {
      if (DOLOG)
	fprintf(logstr, "  MSOLVE: call mmodify and return\n");
      mmodify();

      /* reset the status vector */
      for (j = 0; j < n; j++)
	if (status[j][0] == 'C')
	  status[j][0] = 'c';
      if (DOLOG) {
	fprintf(logstr, "           status=");
        for (j = 0; j < n; j++)
          fprintf(logstr, "%3.3s", status[j]);
	fprintf(logstr, "\n");
      }
      return;
    }

    if (DOLOG)
      fprintf(logstr, "  MSOLVE: isolated %d roots excep=%d\n", nzc, excep);

    if (excep) {
      for (i = 0; i <= n; i++)
	oldpunt[i] = punt[i];
      oldnclust = nclust;

      /* cluster analysis */
      if (DOLOG)
	fprintf(logstr, "  MSOLVE: call mcluster\n");

      mcluster(2 * n); /* Isolation factor */
      
      newtis_old = newtis;
      if(newtis == 0) 
	mnewtis();
      if(DOLOG)
      	fprintf(logstr,"  MSOLVE: newtis_old=%d, newtis=%d, oldncl=%d, nclust=%d\n",
      	                newtis_old, newtis, oldnclust, nclust);

      if (oldnclust == nclust && !(newtis == 1 && newtis_old == 0))
/*#		if(&& iter != 0) AGO99 */
	{
/*#D !newtis */
	  if (DOLOG)
	    fprintf(logstr, "  MSOLVE: CYCLE\n");
	  continue;
	} else {
	  if (DOLOG)
	    fprintf(logstr, "  MSOLVE: call modify\n");
	  mmodify();
	  if (iter == 0)
	    /* if first packet: reset the status vector */
	    for (j = 0; j < n; j++)
	      if (status[j][0] == 'C')
		status[j][0] = 'c';
	  if (DOLOG) {
	    fprintf(logstr, "  MSOLVE:  status=");
            for (j = 0; j < n; j++)
              fprintf(logstr, "%3.3s", status[j]);
            fprintf(logstr, "\n");
          }
	  /* If the polynomial is not given in terms of its coeff. then
	   * skip the restart stage */
	  if (data_type[0] != 'u') {
	  /* choose new starting approximations only for new clusters */
	  if (DOLOG)
	    fprintf(logstr, "  MSOLVE: call mrestart for new clusters\n");
	  mrestart();
	}
	/* reset the status vector */
	for (j = 0; j < n; j++) {
	  if (status[j][0] == 'C')
	    status[j][0] = 'c';
	  again_old[j] = again[j];
	}
	/* update 'again' */
	if (DOLOG)
	  fprintf(logstr, "  MSOLVE: call update2 : ");
	update();
	if (DOLOG) {
	  fprintf(logstr, "  MSOLVE: again = ");
	  for (j = 0; j < n; j++)
	    fprintf(logstr, "%d", again[j]);
	  fprintf(logstr, "\n");
	}
	/* adjust 'again'  This is needed since we are between two packets */
	for (i = 0; i < n; i++)
	  if (!again_old[i])
	    again[i] = false;
	if (DOLOG) {
	  fprintf(logstr, "  MSOLVE: adjusted again = ");
	  for (j = 0; j < n; j++)
	    fprintf(logstr, "%d", again[j]);
	  fprintf(logstr, "\n");
	}
	if (DOLOG)
	  fprintf(logstr, "  MSOLVE: call checkstop\n");
	if (check_stop()) {
	  mmodify();		

	  /* reset the status vector */
	  for (j = 0; j < n; j++)
	    if (status[j][0] == 'C')
	      status[j][0] = 'c';

	  return;
	}

	nzc = 0;
	if (goal[1] == 'a' && goal[2] == 'n' && goal[3] == 'n')
	  for (i = 0; i < n; i++)
	    if (status[i][0] == 'i' || status[i][0] == 'a' || status[i][0] == 'o')
	      nzc++;
	if (DOLOG)
	  fprintf(logstr, "  MSOLVE: check again nzc=%d\n", nzc);
	if (nzc == n) {
	  if (DOLOG)
	    fprintf(logstr, "  MSOLVE: call mmodify and return");
	  mmodify();

	  /* reset the status vector */
	  for (j = 0; j < n; j++)
	    if (status[j][0] == 'C')
	      status[j][0] = 'c';

	  return;
	}
      }
    } else
      break;
  }

  if (iter == max_pack) {
    error(1, "MP: reached the maximum number of packet iteration");
  }
  
  if (DOLOG) {
    fprintf(logstr, "  MSOLVE: MP: nit= %d\n", nit);
    fprintf(logstr, "  MSOLVE: call mcluster\n");
  }
  
      mcluster(2 * n); /* Isolation factor */

  if (DOLOG)
    fprintf(logstr, "  MSOLVE:  call mmodify\n");
  mmodify();

  for (j = 0; j < n; j++)
    if (status[j][0] == 'C')
      status[j][0] = 'c';

  if (DOLOG) {
    fprintf(logstr, "  MSOLVE: status=");
    for (j = 0; j < n; j++)
      fprintf(logstr, "%3.3s", status[j]);
    fprintf(logstr, "\n");
    fprintf(logstr, "  MSOLVE: call update3 : ");
  }
  update();
  
  if (DOLOG) {
    for (j = 0; j < n; j++)
      fprintf(logstr, "%d", again[j]);
    fprintf(logstr, "\n");
  }
}

/*********************************************************
*                SUBROUTINE MPOLZER                      *
**********************************************************
 This routine applies nit iterations of Aberth's method
 to the i-th component of the approximations for which
 again[i] is true
 Set again[i]=false if the i-th approximation is in the
 root neighbourhood.  Stop if  again[i]=false for any i.
 excep= true if after nit iterations some approximation
 is still out of the root neighbourhood.
 ********************************************************/
void
mpolzer(int *it, boolean * excep)
{
  int nzeros, i, j, iter, l;
  tmpc_t corr, abcorr;
  rdpe_t eps, rad1, rtmp;
  cdpe_t ctmp;

  tmpc_init2(abcorr, mpwp);
  tmpc_init2(corr, mpwp);

  rdpe_mul_d(eps, mp_epsilon, (double) 4 * n);

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  nzeros = 0;
  for (i = 0; i < n; i++)
    if (!again[i])
      nzeros++;
  if (nzeros == n)
    goto endfun;

  /* Start Aberth's iterations */
  for (iter = 0; iter < max_it; iter++) {	/* do_iter: */
    for (j = 0; j < nclust; j++) {	/* do_clust: */
      for (i = 0; i < punt[j + 1] - punt[j]; i++) {	/* do_indice: */
	l = clust[punt[j] + i];
	if (again[l]) {
	  (*it)++;
	  if (data_type[0] != 'u') {
	    /* sparse/dense polynomial */
	    rdpe_set(rad1, drad[l]);
	    mnewton(n, mroot[l], drad[l], corr, mfpc, mfppc,
		    dap, spar, &again[l]);
	    if (iter == 0 && !again[l] && rdpe_gt(drad[l], rad1) &&
		rdpe_ne(rad1, rdpe_zero))
	      rdpe_set(drad[l], rad1);

      /************************************************
        The above condition is needed to cope with the case
        where at the first iteration the starting point is
        already in the root neighbourhood and the actually
        computed radius is too big since the value of the
        first derivative is too small.
        In this case the previous radius bound, obtained by
        means of Rouche' is more reliable and strict
        ***********************************************/
	  } else		/* user's polynomial */
	    mnewton_usr(mroot[l], drad[l], corr, &again[l]);

	  if (again[l] ||
	      /* the correction is performed only if iter!=1 or rad[l]!=rad1 */
	      data_type[0] == 'u' || iter != 0 || rdpe_ne(drad[l], rad1)) {
	    maberth_s(l, j, abcorr);
	    mpc_mul_eq(abcorr, corr);
	    mpc_neg_eq(abcorr);
	    mpc_add_eq_ui(abcorr, 1, 0);
	    mpc_div(abcorr, corr, abcorr);
	    mpc_sub_eq(mroot[l], abcorr);
	    mpc_get_cdpe(ctmp, abcorr);
	    cdpe_mod(rtmp, ctmp);
	    rdpe_add_eq(drad[l], rtmp);
	  }

	/* check for new approximated roots */
	  if (!again[l]) {
	    nzeros++;
	    if (nzeros == n)
	      goto endfun;
	  }

	}
      }
    }
    if (nzeros == n)
      goto endfun;
  }
  *excep = true;

endfun:			/* free local MP variables */
  tmpc_clear(corr);
  tmpc_clear(abcorr);
}
