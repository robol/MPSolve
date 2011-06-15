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

/********************************************************
*      SUBROUTINE INCLUSION                             *
*********************************************************/
mps_boolean
mps_inclusion(mps_status* s)
{
  int i, j, k, n1, n2, oldnclust;
  rdpe_t rad, difr;
  cdpe_t difc;
  tmpc_t tmp;


  rdpe_t ap, az, temp, ep, apeps;
  cdpe_t temp1;
  tmpc_t p;
  

  /* add inclusion code here */
  if (!s->chkrad || s->lastphase != mp_phase) {
    if (s->DOLOG)
      fprintf(s->logstr, "Skipping inclusion disks check.\n");
    return true;
  }

  if (s->DOLOG)
    fprintf(s->logstr, "Checking inclusion disks...\n");

  if (s->DOLOG) {
    fprintf(s->logstr, "Old radii\n");
    for (i=0; i<s->n; i++) {
      fprintf(s->logstr, "r(%d)=", i);
      rdpe_outln_str(s->logstr, s->drad[i]);
    }
  }

  /* save old radii */
  for (i=0; i<s->n; i++)
    rdpe_set(s->dap1[i], s->drad[i]);

  tmpc_init2(p, s->mpwp);
  rdpe_mul_d(ep, s->mp_epsilon, (double) (s->n * 4));
  
  tmpc_init2(tmp, s->mpwp);

  for (i=0; i<s->n; i++) {
    
    /* compute denominator */
    rdpe_set(rad, rdpe_one);
    for (j=0; j<s->n; j++) {
      if (i==j) continue;
      mpc_sub(tmp, s->mroot[j], s->mroot[i]);
      mpc_get_cdpe(difc, tmp);
      cdpe_smod(difr, difc);
      rdpe_mul_eq(rad, difr);
    }
    rdpe_sqrt_eq(rad);
    rdpe_mul_eq(rad, s->dap[s->n]);
    
    /* compute numerator*/
    if (s->data_type[0] == 's') {	/* case of sparse polynomial */

      n1 = s->n + 1;
      n2 = s->n;
      
      /* compute p(mroot[i]) */
      mps_parhorner(s, n1, s->mroot[i], s->mfpc, s->spar, p);
      mpc_get_cdpe(temp1, s->mroot[i]);
      cdpe_mod(az, temp1);
      
      /* compute bound to the error */
      mps_aparhorner(s, n1, az, s->dap, s->spar, ap);
      
    } else {			/*  dense polynomial */
      
      /* commpute p(mroot[i]) and p'(mroot[i]) */
      mpc_set(p, s->mfpc[s->n]);
      for (k = s->n - 1; k > 0; k--) {
	mpc_mul(p, p, s->mroot[i]);
	mpc_add(p, p, s->mfpc[k]);
      }
      mpc_mul(p, p, s->mroot[i]);
      mpc_add(p, p, s->mfpc[0]);
      
      /* compute bound to the error */
      rdpe_set(ap, s->dap[s->n]);
      mpc_get_cdpe(temp1, s->mroot[i]);
      cdpe_mod(az, temp1);
      for (k = s->n - 1; k >= 0; k--) {
	rdpe_mul(temp, ap, az);
	rdpe_add(ap, temp, s->dap[k]);
      }
    }
    
    /* common part */
    mpc_get_cdpe(difc, p);
    cdpe_mod(difr, difc);
    rdpe_mul(apeps, ap, ep);
    rdpe_add_eq(apeps, difr);
    rdpe_mul_eq_d(apeps, (double) s->n);

    /* compute ratio */
    rdpe_div(s->drad[i], apeps, rad);
    
    if (s->DOLOG) {
      fprintf(s->logstr, "New r(%d)=", i);
      rdpe_outln_str(s->logstr, s->drad[i]);
    }
  }
 
  oldnclust = s->nclust;

  mps_mcluster(s, 2*s->n);

  if (s->nclust>=oldnclust) {
    /* choose the smallest radius */
    for (i=0; i<s->n; i++)
      if (rdpe_lt(s->dap1[i], s->drad[i]))
	rdpe_set(s->drad[i], s->dap1[i]);
    /* update(); */
  } else
    mps_warn(s, "Some roots might be not approximated");

  tmpc_clear(tmp);  
  tmpc_clear(p);

  return true;
}
