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

#include <mps/mps_poly.h>
#include <mps/rursolve.h>

mpz_t *mpdemo = NULL;	/* imaginary part of the integer input coeff. */

/***********************************************************
 *           SUBROUTINE RURSOLVE                           *
 ***********************************************************
  reads and computes the univariate rational representation
 ***********************************************************/
void
rursolve(mps_status* s)
{
  mpspoly_t p;
  mpz_t tnden;
  mpf_t rn, rd, mptemp;
  mpc_t crn, crd, res;
  mpc_t **sol;
  int **precsol;
  int i, j, pos, m, degmax, dig, dig2, deg, iprec, dprec, nc, tn, read_bytes;

  /* ==============================
   * phase 1: input global values
   * ============================== */

  /* Assign to i to get around warnings */
  read_bytes = fscanf(s->instr, "%d", &m);
  read_bytes = fscanf(s->instr, "%d", &degmax);
  read_bytes = fscanf(s->instr, "%d", &dig);
  dig2 = (int) ((dig + 1) * LOG2_10);

  /* ==============================
   * phase 2: read and solve univariate problem
   * ============================== */

  /* read polynomial */
  mps_read_poly(s, s->instr, p);
  
  /* set polynomial data */
  mps_set_poly(s, p);

  /* allocate global variables */
  mps_allocate_data(s);

  if (s->DOLOG)
    fprintf(s->logstr, "Calling mpsolve...\n");

  mps_mpsolve(s);
  
  /* copy computed roots */
  mps_copy_roots(s);

  /* ==============================
   * phase 3: allocate rur variables
   * ============================== */

  mpz_init(tnden);
  mpf_init2(rn, dig2);
  mpf_init2(rd, dig2);
  mpf_init2(mptemp, dig2);
  mpc_init2(crn, dig2);
  mpc_init2(crd, dig2);
  mpc_init2(res, dig2);

  /* allocate matrices */
  tn = (s->zero_roots) ? s->n + 1 : s->n;

  sol = (mpc_t **) malloc(tn * sizeof(mpc_t *));
  for (i = 0; i < tn; i++)
    sol[i] = (mpc_t *) malloc(m * sizeof(mpc_t));
  for (i = 0; i < tn; i++)
    for (j = 0; j < m; j++)
      mpc_init2(sol[i][j], dig2);

  precsol = (int **) malloc(tn * sizeof(int *));
  for (i = 0; i < tn; i++)
    precsol[i] = (int *) malloc(m * sizeof(int));

  /* mpdemo */
  mpdemo = mpz_valloc(degmax + 1);
  for (i = 0; i <= degmax; i++)
    mpz_init(mpdemo[i]);

  /* recycle mpsolve variables */
  if (degmax > s->n) {
    s->dap1 = (rdpe_t *) realloc(s->dap1, (degmax + 1) * sizeof(rdpe_t));
    s->mfpc1 = (mpc_t *) realloc(s->mfpc1, (degmax + 1) * sizeof(mpc_t));
    s->mfpc2 = (mpc_t *) realloc(s->mfpc2, (degmax + 1) * sizeof(mpc_t));
    for (i = s->n; i <= degmax; i++) {
      rdpe_set(s->dap1[i], rdpe_zero);
      mpc_init(s->mfpc1[i]);
      mpc_init(s->mfpc2[i]);
    }
  }
  
  /* ==============================
   * Phase 4: Read denominator polynomial
   * ============================== */

  read_bytes = fscanf(s->instr, "%3s", s->data_type);
  read_bytes = fscanf(s->instr, "%d", &deg);

  if (s->data_type[0] == 'd') {	/* dense polynomial */
    for (i = 0; i <= deg; i++) {
      mpz_inp_str(mpdemo[i], s->instr, 10);
      /* represent integers as dpe */
      mpf_set_z(mptemp, mpdemo[i]);
      mpf_get_rdpe(s->dap1[i], mptemp);
      rdpe_abs_eq(s->dap1[i]);
    }
  } else {			/* sparse polynomial  */
    for (i = 0; i <= deg; i++) {
      rdpe_set(s->dap1[i], rdpe_zero);
      mpz_set_ui(mpdemo[i], 0);
    }
    read_bytes = fscanf(s->instr, "%d", &nc);
    for (i = 0; i < nc; i++) {
      read_bytes = fscanf(s->instr, "%d", &pos);
      mpz_inp_str(mpdemo[pos], s->instr, 10);
      /* represent integers as dpe */
      mpf_set_z(mptemp, mpdemo[pos]);
      mpf_get_rdpe(s->dap1[pos], mptemp);
      rdpe_abs_eq(s->dap1[pos]);
    }
  }

  mpz_set(tnden, mpdemo[0]);

  /* ==============================
   * Phase 5: compute denominators
   * ============================== */

  for (i = 0; i < s->n; i++) {
    mps_horner(s, res, &dprec, &iprec, deg, i);

    /* Computed:
     * iprec: bit precision of mroot[i]
     * dprec: bit precision of the denominator mfpc[i] */

    if (s->DOLOG) {
      fprintf(s->logstr, " Computed denominator=");
      mpc_out_str(s->logstr, 10, 0, res);
      fprintf(s->logstr, "\n i=%d, Prec root=%d, prec denom=%d\n", i, iprec, dprec);
    }

    /* check precision */
    if (dprec < dig2) {
      mps_refine(s, i, iprec + dig2 - dprec);
      mps_horner(s, res, &dprec, &iprec, deg, i);

      if (s->DOLOG) {
	fprintf(s->logstr, " Computed denominator=");
	mpc_out_str(s->logstr, 10, 0, res);
	fprintf(s->logstr, "\n i=%d, Prec root=%d, prec denom=%d\n", i, iprec, dprec);
      }
      if (s->DOLOG) {
	fprintf(s->logstr, " Refined denominatorn=");
	mpc_out_str(s->logstr, 10, 0, res);
	fprintf(s->logstr, "\n");
      }
    }
    
    /* store denominator */
    mpc_set(sol[i][m - 1], res);
    precsol[i][m - 1] = dprec;

    if (s->DOLOG) {
      fprintf(s->logstr, "Denominator[%d] = ", i);
      mpc_out_str(s->logstr, 10, 0, sol[i][m - 1]);
      fprintf(s->logstr, "\n");
    }
  }

  /* ==============================
   * Phase 6: compute all numerators
   * ============================== */

  for (j = 0; j < m; j++) {	/* main loop */
    /* Read j-th numerator polynomial */
    read_bytes = fscanf(s->instr, "%3s", s->data_type);
    read_bytes = fscanf(s->instr, "%d", &deg);

    if (s->data_type[0] == 'd') {	/* dense polynomial */
      for (i = 0; i <= deg; i++) {
	mpz_inp_str(mpdemo[i], s->instr, 10);
	/* represent integers as dpe */
	mpf_set_z(mptemp, mpdemo[i]);
	mpf_get_rdpe(s->dap1[i], mptemp);
	rdpe_abs_eq(s->dap1[i]);
      }
    } else {			/* sparse polynomial */
      for (i = 0; i <= deg; i++) {
	rdpe_set(s->dap1[i], rdpe_zero);
	mpz_set_ui(mpdemo[i], 0);
      }
      read_bytes = fscanf(s->instr, "%d", &nc);
      for (i = 0; i < nc; i++) {
	read_bytes = fscanf(s->instr, "%d", &pos);
	mpz_inp_str(mpdemo[pos], s->instr, 10);
	/* represent integers as dpe */
	mpf_set_z(mptemp, mpdemo[pos]);
	mpf_get_rdpe(s->dap1[pos], mptemp);
	rdpe_abs_eq(s->dap1[pos]);
      }
    }

    /* compute numerators */
    
    for (i = 0; i < s->n; i++) {
      mps_horner(s, res, &dprec, &iprec, deg, i);

      if (s->DOLOG)
	fprintf(s->logstr, "i=%d, Prec zero=%d, prec num=%d\n", i, iprec, dprec);

      if (dprec < dig2) {
	mps_refine(s, i, iprec + dig2 - dprec);
	mps_horner(s, res, &dprec, &iprec, deg, i);
	if (s->DOLOG)
	  fprintf(s->logstr, "i=%d, Prec zero=%d, prec num=%d\n", i, iprec, dprec);
      }
      
      if (s->DOLOG) {
	fprintf(s->logstr, " Numerator (%d, %d) = ", j, i);
	mpc_out_str(s->logstr, 10, 0, res);
	fprintf(s->logstr, "\n");
      }

      /* compute and store root i, j*/
      mpc_div(sol[i][j], res, sol[i][m - 1]);
      precsol[i][j] = MIN(dprec, precsol[i][m - 1]);

      if (s->DOLOG)
	fprintf(s->logstr, "Root prec. i=%d, j=%d prec=%d\n", i, j, precsol[i][j]);
    }

    if (s->zero_roots) {
      mpf_set_z(rn, mpdemo[0]);
      mpf_set_z(rd, tnden);
      mpf_div(rn, rn, rd);
      mpf_set(mpc_Re(crn), rn);
      mpf_set_ui(mpc_Im(crn), 0);
      mpc_set(sol[s->n][j], crn);
    }
  }				/* main loop */

  /* ==============================
   * Phase 7: output
   * ============================== */

  for (pos = 0; pos < s->n; pos++) {
    i = s->order[pos];
    fprintf(s->outstr, "\nRoot(%d) = (", pos + 1);
    for (j = 0; j < m; j++) {
      mps_ruroutroot(s, sol[i][j], s->status[i][1], precsol[i][j], dig2);
      if (j < m - 1)
	fprintf(s->outstr, ",  ");
    }
    fprintf(s->outstr, ")\n");
  }

  if (s->n < tn) {
    fprintf(s->outstr, "\nRoot(%d) = (", s->n + 1);
    for (j = 0; j < m; j++) {
      mpf_set(rn, mpc_Re(sol[s->n][j]));
      mpf_out_str(s->outstr, 10, dig, rn);
      if (j < m - 1)
	fprintf(s->outstr, ",  ");
    }
    fprintf(s->outstr, ")\n");
  }
  
  /* ==============================
   * Phase 8: free dynamic variables
   * ============================== */

  for (i = 0; i < tn; i++)
    for (j = 0; j < m; j++)
      mpc_clear(sol[i][j]);
  for (i = 0; i < tn; i++)
    free(sol[i]);
  free(sol);

  for (i = 0; i < tn; i++)
    free(precsol[i]);
  free(precsol);

  mpc_clear(res);
  mpc_clear(crn);
  mpc_clear(crd);
  mpf_clear(rn);
  mpf_clear(rd);
  mpf_clear(mptemp);
  mpz_clear(tnden);
}

/*********************************************************
*      SUBROUTINE RUROUTROOT                           *
*********************************************************/
void
mps_ruroutroot(mps_status* s, mpc_t root, char status, long prec, long out_prec)
{
  tmpf_t t;
  cdpe_t ctmp;
  rdpe_t rtmp;
  long l, lr, outdig;

  prec = MIN(prec, out_prec);
  tmpf_init2(t, prec);

  mpc_get_cdpe(ctmp, root);
  cdpe_mod(rtmp, ctmp);
  l = rdpe_Esp(rtmp);

  mpf_get_rdpe(rtmp, mpc_Re(root));
  lr = rdpe_Esp(rtmp);
  if (lr == 0 && rdpe_Mnt(rtmp) == 0)
    fprintf(s->outstr, "0e%ld", (long) (-prec * LOG10_2));
  else {
    outdig = prec - (l - lr);
    if (outdig <= 0)
      fprintf(s->outstr, "0.e%ld", (long) (lr * LOG10_2));
    else {
      outdig = (long) (outdig * LOG10_2);
      if (mpf_sgn(mpc_Re(root)) < 0)
	fprintf(s->outstr, "-");
      mpf_abs(t, mpc_Re(root));
      mpf_out_str(s->outstr, 10, outdig, t);
    }
  }

  if (status != 'R') {
    if (mpf_sgn(mpc_Im(root)) >= 0)
      fprintf(s->outstr, " + I * ");
    else
      fprintf(s->outstr, " - I * ");
    mpf_get_rdpe(rtmp, mpc_Im(root));
    lr = rdpe_Esp(rtmp);
    outdig = prec - (l - lr);
    if (outdig <= 0)
      fprintf(s->outstr, "0.e%ld", (long) (lr * LOG10_2));
    else {
      outdig = (long) (outdig * LOG10_2);
      mpf_abs(t, mpc_Im(root));
      mpf_out_str(s->outstr, 10, outdig, t);
    }
  }
  tmpf_clear(t);
}
