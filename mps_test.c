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

/********************************************************
*      SUBROUTINE INCLUSION                             *
*********************************************************/
boolean
inclusion(void)
{
  int i, j, k, n1, n2, oldnclust;
  rdpe_t rad, difr;
  cdpe_t difc;
  tmpc_t tmp;


  rdpe_t ap, az, temp, ep, apeps;
  cdpe_t temp1;
  tmpc_t p;
  

  /* add inclusion code here */
  if (!chkrad || lastphase != mp_phase) {
    if (DOLOG)
      fprintf(logstr, "Skipping inclusion disks check.\n");
    return true;
  }

  if (DOLOG)
    fprintf(logstr, "Checking inclusion disks...\n");

  if (DOLOG) {
    fprintf(logstr, "Old radii\n");
    for (i=0; i<n; i++) {
      fprintf(logstr, "r(%d)=", i);
      rdpe_outln_str(logstr, drad[i]);
    }
  }

  /* save old radii */
  for (i=0; i<n; i++)
    rdpe_set(dap1[i], drad[i]);

  tmpc_init2(p, mpwp);
  rdpe_mul_d(ep, mp_epsilon, (double) (n * 4));
  
  tmpc_init2(tmp, mpwp);

  for (i=0; i<n; i++) {
    
    /* compute denominator */
    rdpe_set(rad, rdpe_one);
    for (j=0; j<n; j++) {
      if (i==j) continue;
      mpc_sub(tmp, mroot[j], mroot[i]);
      mpc_get_cdpe(difc, tmp);
      cdpe_smod(difr, difc);
      rdpe_mul_eq(rad, difr);
    }
    rdpe_sqrt_eq(rad);
    rdpe_mul_eq(rad, dap[n]);
    
    /* compute numerator*/
    if (data_type[0] == 's') {	/* case of sparse polynomial */

      n1 = n + 1;
      n2 = n;
      
      /* compute p(mroot[i]) */
      parhorner(n1, mroot[i], mfpc, spar, p);
      mpc_get_cdpe(temp1, mroot[i]);
      cdpe_mod(az, temp1);
      
      /* compute bound to the error */
      aparhorner(n1, az, dap, spar, ap);
      
    } else {			/*  dense polynomial */
      
      /* commpute p(mroot[i]) and p'(mroot[i]) */
      mpc_set(p, mfpc[n]);
      for (k = n - 1; k > 0; k--) {
	mpc_mul(p, p, mroot[i]);
	mpc_add(p, p, mfpc[k]);
      }
      mpc_mul(p, p, mroot[i]);
      mpc_add(p, p, mfpc[0]);
      
      /* compute bound to the error */
      rdpe_set(ap, dap[n]);
      mpc_get_cdpe(temp1, mroot[i]);
      cdpe_mod(az, temp1);
      for (k = n - 1; k >= 0; k--) {
	rdpe_mul(temp, ap, az);
	rdpe_add(ap, temp, dap[k]);
      }
    }
    
    /* common part */
    mpc_get_cdpe(difc, p);
    cdpe_mod(difr, difc);
    rdpe_mul(apeps, ap, ep);
    rdpe_add_eq(apeps, difr);
    rdpe_mul_eq_d(apeps, (double) n);

    /* compute ratio */
    rdpe_div(drad[i], apeps, rad);
    
    if (DOLOG) {
      fprintf(logstr, "New r(%d)=", i);
      rdpe_outln_str(logstr, drad[i]);
    }
  }
 
  oldnclust = nclust;

  mcluster(2*n);

  if (nclust>=oldnclust) {
    /* choose the smallest radius */
    for (i=0; i<n; i++)
      if (rdpe_lt(dap1[i], drad[i]))
	rdpe_set(drad[i], dap1[i]);
    /* update(); */
  } else
    warn("Some roots might be not approximated");

  tmpc_clear(tmp);  
  tmpc_clear(p);

  return true;
}
