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

#include <stdarg.h>
#include <string.h>
#include <mps/core.h>

#define ISZERO -1

/*********************************************************
*      SUBROUTINE READROOTS                              *
*********************************************************/
void
msp_readroots(mps_status* s)
{
  long digits;
  int i, read_elements;

  if (s->DOLOG)
    fprintf(s->logstr, "Reading roots...\n");

  read_elements = fscanf(s->rtstr, "%ld", &digits);
  if (!read_elements) {
    mps_error(s, 1, "Error while reading roots, aborting.");
  }

  /* precision setup code goes here */

  for(i=0; i<s->n; i++)
    mpc_inp_str_u(s->mroot[i], s->rtstr, 10);
}

/*********************************************************
*      SUBROUTINE COUNTROOTS                             *
*********************************************************/
void
mps_countroots(mps_status* s)
{
  int k;

  if (s->DOLOG)
    fprintf(s->logstr, "Counting roots...\n");

  s->count[0] = s->count[1] = s->count[2] = 0;

  for (k = 0; k < s->n; k++)
    switch (s->status[k][2]) {
    case 'i':
      s->count[0]++;
      break;
    case 'o':
      s->count[1]++;
      break;
    default:
      s->count[2]++;
      break;
    }

  if (s->goal[1] == 'o')
    s->count[1] += s->zero_roots;
  else
    s->count[0] += s->zero_roots;
}

/*********************************************************
*      SUBROUTINE OUTCOUNT                               *
*********************************************************/
void
mps_outcount(mps_status* s)
{
  mps_countroots(s);

  fprintf(s->outstr, "%d roots are inside;\n", s->count[0]);
  fprintf(s->outstr, "%d roots are outside;\n", s->count[1]);
  fprintf(s->outstr, "%d roots are uncertain.\n", s->count[2]);
  if (s->DOLOG) {
    fprintf(s->logstr, "%d roots are inside;\n", s->count[0]);
    fprintf(s->logstr, "%d roots are outside;\n", s->count[1]);
    fprintf(s->logstr, "%d roots are uncertain.\n", s->count[2]);
  }
}

/*********************************************************
*      SUBROUTINE OUTFLOAT                               *
*********************************************************/
void
mps_outfloat(mps_status* s, mpf_t f, rdpe_t rad, long out_digit, mps_boolean sign)
{
  tmpf_t t;
  rdpe_t r, ro;
  double d;
  long l, digit, true_digit;

  if (s->goal[4] == 'f') {
    tmpf_init2(t, mpf_get_prec(f));
    mpf_set(t, f);
    mpf_out_str(s->outstr, 10, 0, t);
    tmpf_clear(t);
    return;
  }

  tmpf_init2(t, s->prec_out);

  mpf_get_rdpe(ro, f);
  if (s->goal[4] == 'g')
    rdpe_out_str_u(s->outstr, ro);
  else {
    rdpe_abs_eq(ro);
    if (rdpe_ne(ro, rdpe_zero))
      rdpe_div(r, rad, ro);
    else
      rdpe_set_d(r, 1.0e-10);
    digit = (long) (-rdpe_log10(r) - 0.5);
    if (digit <= 0) {
      rdpe_get_dl(&d, &l, ro);
      fprintf(s->outstr, "0.e%ld", l);
    } else {
      true_digit = (long) (LOG10_2 * mpf_get_prec(f));
      true_digit = MIN(digit, true_digit);
      true_digit = MIN(true_digit, out_digit);
      if (sign)
	mpf_set(t, f);
      else
	mpf_abs(t, f);
      mpf_out_str(s->outstr, 10, true_digit, t);
    }
  }

  tmpf_clear(t);
}

/*********************************************************
*      SUBROUTINE OUTROOT                                *
*********************************************************/
void
mps_outroot(mps_status* s, int i)
{
  static int num = 0; /* output roots count */
  long out_digit;

  out_digit = (long) (LOG10_2 * s->prec_out) + 10;
  num++;

  /* print format header */
  switch (s->goal[4]) {
  case 'c':
  case 'f':
    fprintf(s->outstr, "(");
    break;
  case 'v':
    fprintf(s->outstr, "Root(%d) = ", num);
    break;
  }

  /* print real part */
  if (i == ISZERO || s->status[i][1] == 'I')
    fprintf(s->outstr, "0");
  else
    mps_outfloat(s, mpc_Re(s->mroot[i]), s->drad[i], out_digit, true);

  /* print format middle part */
  switch (s->goal[4]) {
  case 'b':
    fprintf(s->outstr, " ");
    break;
  case 'g':
    fprintf(s->outstr, "\t");
    break;
  case 'c':
  case 'f':
    fprintf(s->outstr, ", ");
    break;
  case 'v':
    if (i == ISZERO || mpf_sgn(mpc_Im(s->mroot[i])) >= 0)
      fprintf(s->outstr, " + I * ");
    else
      fprintf(s->outstr, " - I * ");
    break;
  }

  /* print imaginary part */
  if (i == ISZERO || s->status[i][1] == 'R')
    fprintf(s->outstr, "0");
  else
    mps_outfloat(s, mpc_Im(s->mroot[i]), s->drad[i], out_digit, s->goal[4] != 'v');

  /* print format ending */
  switch (s->goal[4]) {
  case 'c':
    fprintf(s->outstr, ")");
    break;
  case 'f':
    fprintf(s->outstr, ")\n");
    if (i != ISZERO) {
      rdpe_outln_str(s->outstr, s->drad[i]);
      fprintf(s->outstr, "%4.3s\n", s->status[i]);
    } else
      fprintf(s->outstr, " 0\n ---\n");
    break;
  }
  fprintf(s->outstr, "\n");

  /* debug info */
  if (s->DOLOG) {
    if(i == ISZERO)
      fprintf(s->logstr, "zero root %-4d = 0", num);
    else {
      fprintf(s->logstr, "root %-4d = ", i);
      mpc_out_str_2(s->logstr, 10, 0, 0, s->mroot[i]);
      fprintf(s->logstr, "\n");
      fprintf(s->logstr, "  radius = ");
      rdpe_outln_str(s->logstr, s->drad[i]);
      fprintf(s->logstr, "  prec = %ld\n",
	      (long) (mpc_get_prec(s->mroot[i])/LOG2_10));
      fprintf(s->logstr, "  status = %4.3s\n", s->status[i]);
      fprintf(s->logstr, "--------------------\n");
    }
  }
}

/*********************************************************
*      SUBROUTINE OUTPUT                                 *
*********************************************************/
void
mps_output(mps_status* s)
{
  int i, ind;

  if (s->DOLOG)
    fprintf(s->logstr, "--------------------\n");

  if (s->goal[0] == 'c')
    mps_outcount(s);
  else {
    if (s->goal[1] != 'o')
      for (i = 0; i < s->zero_roots; i++)
        mps_outroot(s, ISZERO);
    for (ind = 0; ind < s->n; ind++) {
      i = s->order[ind];
      if (s->status[i][2] == 'o')
	continue;
      mps_outroot(s, i);
    }
  }
}

/*********************************************************
*      SUBROUTINE COPY_ROOTS                             *
*********************************************************/
void
mps_copy_roots(mps_status* s)
{
  int i;

  switch (s->lastphase) {
  case no_phase:
    mps_error(s, 1, "Nothing to copy");
    break;

  case float_phase:
    if (s->DOSORT)
      mps_fsort(s);
    for (i = 0; i < s->n; i++) {
      mpc_set_prec(s->mroot[i], DBL_MANT_DIG);
      mpc_set_cplx(s->mroot[i], s->froot[i]);
      rdpe_set_d(s->drad[i], s->frad[i]);
    }
    break;

  case dpe_phase:
    if (s->DOSORT)
      mps_dsort(s);
    for (i = 0; i < s->n; i++) {
      mpc_set_prec(s->mroot[i], DBL_MANT_DIG);
      mpc_set_cdpe(s->mroot[i], s->droot[i]);
    }
    break;

  case mp_phase:
    if (s->DOSORT)
      mps_msort(s);
    break;

  }
}

/*************************************************************
 *                     SUBROUTINE DUMP                       *
 *************************************************************/
void
mps_dump(mps_status* s, FILE * dmpstr)
{
  int i;

  fprintf(dmpstr, "\nDumping...\n");

  /* output current status */
  fprintf(dmpstr,
	  "Phase=%d, In=%d, Out=%d, Uncertain=%d, Zero=%d, Clusters=%d\n",
	  s->lastphase, s->count[0], s->count[1], s->count[2], s->zero_roots, s->nclust);

  /* output current approximations */
  fprintf(dmpstr, "\nCurrent approximations:\n");
  for (i = 0; i < s->n; i++) {
    fprintf(dmpstr, "%d:\t", i);

    switch (s->lastphase) {
    case no_phase:
    case float_phase:
      cplx_outln_str(dmpstr, s->froot[i]);
      break;

    case dpe_phase:
      cdpe_outln_str(dmpstr, s->droot[i]);
      break;

    case mp_phase:
      mpc_outln_str(dmpstr, 10, 0, s->mroot[i]);
      break;
    }
  }

  /* output radii */
  fprintf(dmpstr, "\nCurrent radii:\n");
  for (i = 0; i < s->n; i++) {
    fprintf(dmpstr, "%d:\t", i);

    switch (s->lastphase) {
    case no_phase:
    case float_phase:
      fprintf(dmpstr, "%e\n", s->frad[i]);
      break;

    case dpe_phase:
    case mp_phase:
      rdpe_outln_str(dmpstr, s->drad[i]);
      break;
    }
  }

  /* output position */
  fprintf(dmpstr, "\nPos:\t");
  for (i = 0; i < s->n; i++)
    fprintf(dmpstr, "%4d", i);

  /* output status information */
  fprintf(dmpstr, "\nStatus:\t");
  for (i = 0; i < s->n; i++)
    fprintf(dmpstr, "%4.3s", s->status[i]);

  /* output cluster information */
  fprintf(dmpstr, "\nClust:\t");
  for (i = 0; i < s->n; i++)
    fprintf(dmpstr, "%4d", s->clust[i]);

  fprintf(dmpstr, "\nPunt:\t");
  for (i = 0; i < s->nclust; i++)
    fprintf(dmpstr, "%4d", s->punt[i]);

  fprintf(dmpstr, "\n\n");
}

/**
 * @brief Dump cluster structure to <code>outstr</code>.
 *
 * @param s the mps_status struct pointer.
 * @param outstr The output stream where the cluster structure
 *  will be dumped.
 */
void
mps_dump_cluster_structure(mps_status* s, FILE* outstr)
{
    int i, j;
    fprintf(outstr, "    MPS_DUMP_CLUSTER_STRUCTURE: Dumping cluster structure\n");

    for(i = 0; i < s->nclust; i++) {
        fprintf(outstr, "     Cluster %d contains %d roots:\n", i, s->punt[i+1] - s->punt[i]);

        /* Dump cluster roots, but not more than 15 for line, to make
         * the output readable. */
        for(j = s->punt[i]; j < s->punt[i+1]; j++) {
            /* Go to a newlint if 15 roots are printed out */
            if ((j - s->punt[i]) % 15 == 0) { fprintf(outstr, "\n       "); }

            printf(" %d", s->clust[j]);
        }

        /* Make space untile the next cluster */
        fprintf(outstr, "\n\n");
    }
}

/**
 * @brief Dump status of all the root approximations
 */
void
mps_dump_status(mps_status* s, FILE* outstr)
{
	int i;
	for(i = 0; i < s->n; i++) {
		fprintf(outstr, "s->status[%d] = ", i);
		fprintf(outstr, "'%c' '%c' '%c'\n", s->status[i][0],
				s->status[i][1], s->status[i][2]);
	}
}

/*************************************************************
 *                     SUBROUTINE WARN                       *
 *************************************************************/
void
mps_warn(mps_status* st, char *s)
{
  if (st->DOWARN) {
    if (s[strlen(s)] == '\n') {
      fprintf(st->logstr, "%s", s);
    }
    else {
    	fprintf(st->logstr, "%s\n", s);
    }
  }
}

/*************************************************************
 *                     SUBROUTINE VAERROR                    *
 *************************************************************/
void
mps_error(mps_status* st, int args, ...)
{
  va_list ap;
  char * s;

  mps_warn(st, "! ");		/* output error message */
  va_start(ap, args);
  while(args--) {
    s = va_arg(ap, char *);
    mps_warn(st, s);		/* output error message */
  }
  va_end(ap);

  mps_dump(st, st->logstr);		/* dump status		*/
  exit(EXIT_FAILURE);	/* exit program         */
}

#if __STDC_VERSION__ < 199901L
#ifndef DISABLE_DEBUG
void
MPS_DEBUG(mps_status* s, const char* templ, ...)
{
    
    va_list ap; 
    if (!s->DOLOG)
        return;
    va_start(ap, templ); 
    gmp_vfprintf(s->logstr, templ, ap); 
    fprintf(s->logstr, "\n"); 
}

void
__MPS_DEBUG(mps_status* s, const char* templ, ...)
{
    va_list ap;
    if (!s->DOLOG)
        return;
    va_start(ap, templ);
    gmp_vfprintf(s->logstr, templ, ap); 
}
#endif
#endif
