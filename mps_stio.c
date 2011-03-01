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
#include "mps.h"

#define ISZERO -1

/*********************************************************
*      SUBROUTINE READROOTS                              *
*********************************************************/
void
readroots(void)
{
  long digits;
  int i;

  if (DOLOG)
    fprintf(logstr, "Reading roots...\n");
  
  fscanf(rtstr, "%ld", &digits);

  /* precision setup code goes here */
  
  for(i=0; i<n; i++)
    mpc_inp_str_u(mroot[i], rtstr, 10);
}

/*********************************************************
*      SUBROUTINE COUNTROOTS                             *
*********************************************************/
void
countroots(void)
{
  int k;

 if (DOLOG)
    fprintf(logstr, "Counting roots...\n");

  count[0] = count[1] = count[2] = 0;

  for (k = 0; k < n; k++)
    switch (status[k][2]) {
    case 'i':
      count[0]++;
      break;
    case 'o':
      count[1]++;
      break;
    default:
      count[2]++;
      break;
    }

  if (goal[1] == 'o')
    count[1] += zero_roots;
  else
    count[0] += zero_roots;
}

/*********************************************************
*      SUBROUTINE OUTCOUNT                               *
*********************************************************/
void
outcount(void)
{
  countroots();
  
  fprintf(outstr, "%d roots are inside;\n", count[0]);
  fprintf(outstr, "%d roots are outside;\n", count[1]);
  fprintf(outstr, "%d roots are uncertain.\n", count[2]);
  if (DOLOG) {
    fprintf(logstr, "%d roots are inside;\n", count[0]);
    fprintf(logstr, "%d roots are outside;\n", count[1]);
    fprintf(logstr, "%d roots are uncertain.\n", count[2]);
  }
}

/*********************************************************
*      SUBROUTINE OUTFLOAT                               *
*********************************************************/
void
outfloat(mpf_t f, rdpe_t rad, long out_digit, boolean sign)
{
  tmpf_t t;
  rdpe_t r, ro;
  double d;
  long l, digit, true_digit;

  if (goal[4] == 'f') {
    tmpf_init2(t, mpf_get_prec(f));
    mpf_set(t, f);
    mpf_out_str(outstr, 10, 0, t);
    tmpf_clear(t);
    return;
  }

  tmpf_init2(t, prec_out);
  
  mpf_get_rdpe(ro, f);
  if (goal[4] == 'g')
    rdpe_out_str_u(outstr, ro);
  else {
    rdpe_abs_eq(ro);
    if (rdpe_ne(ro, rdpe_zero))
      rdpe_div(r, rad, ro);
    else
      rdpe_set_d(r, 1.0e-10);
    digit = (long) (-rdpe_log10(r) - 0.5);
    if (digit <= 0) {
      rdpe_get_dl(&d, &l, ro);
      fprintf(outstr, "0.e%ld", l);
    } else {
      true_digit = (long) (LOG10_2 * mpf_get_prec(f));
      true_digit = MIN(digit, true_digit);
      true_digit = MIN(true_digit, out_digit);
      if (sign)
	mpf_set(t, f);
      else
	mpf_abs(t, f);
      mpf_out_str(outstr, 10, true_digit, t);
    }
  }

  tmpf_clear(t);
}

/*********************************************************
*      SUBROUTINE OUTROOT                                *
*********************************************************/
void
outroot(int i)
{
  static int num = 0; /* output roots count */
  long out_digit;
  
  out_digit = (long) (LOG10_2 * prec_out) + 10;
  num++;

  /* print format header */
  switch (goal[4]) {
  case 'c':
  case 'f':
    fprintf(outstr, "(");
    break;
  case 'v':
    fprintf(outstr, "Root(%d) = ", num);
    break;
  }

  /* print real part */
  if (i == ISZERO || status[i][1] == 'I')
    fprintf(outstr, "0");
  else
    outfloat(mpc_Re(mroot[i]), drad[i], out_digit, true);

  /* print format middle part */
  switch (goal[4]) {
  case 'b':
    fprintf(outstr, " ");
    break;
  case 'g':
    fprintf(outstr, "\t");
    break;
  case 'c':
  case 'f':
    fprintf(outstr, ", ");
    break;
  case 'v':
    if (i == ISZERO || mpf_sgn(mpc_Im(mroot[i])) >= 0)
      fprintf(outstr, " + I * ");
    else
      fprintf(outstr, " - I * ");
    break;
  }

  /* print imaginary part */
  if (i == ISZERO || status[i][1] == 'R')
    fprintf(outstr, "0");
  else
    outfloat(mpc_Im(mroot[i]), drad[i], out_digit, goal[4] != 'v');

  /* print format ending */
  switch (goal[4]) {
  case 'c':
    fprintf(outstr, ")");
    break;
  case 'f':
    fprintf(outstr, ")\n");
    if (i != ISZERO) {
      rdpe_outln_str(outstr, drad[i]);
      fprintf(outstr, "%4.3s\n", status[i]);
    } else
      fprintf(outstr, " 0\n ---\n");
    break;
  }
  fprintf(outstr, "\n");

  /* debug info */
  if (DOLOG) {
    if(i == ISZERO)
      fprintf(logstr, "zero root %-4d = 0", num);
    else {
      fprintf(logstr, "root %-4d = ", i);
      mpc_out_str_2(logstr, 10, 0, 0, mroot[i]);
      fprintf(logstr, "\n");
      fprintf(logstr, "  radius = ");
      rdpe_outln_str(logstr, drad[i]);
      fprintf(logstr, "  prec = %ld\n", 
	      (long) (mpc_get_prec(mroot[i])/LOG2_10));
      fprintf(logstr, "  status = %4.3s\n", status[i]);
      fprintf(logstr, "--------------------\n");
    }
  }
}

/*********************************************************
*      SUBROUTINE OUTPUT                                 *
*********************************************************/
void
output(void)
{
  int i, ind;

  if (DOLOG)
    fprintf(logstr, "--------------------\n");

  if (goal[0] == 'c')
    outcount();
  else {
    if (goal[1] != 'o')
      for (i = 0; i < zero_roots; i++)
        outroot(ISZERO);
    for (ind = 0; ind < n; ind++) {
      i = order[ind];
      if (status[i][2] == 'o')
	continue;
      outroot(i);
    }
  }
}

/*********************************************************
*      SUBROUTINE COPY_ROOTS                             *
*********************************************************/
void
copy_roots(void)
{
  int i;

  switch (lastphase) {
  case no_phase:
    error(1, "Nothing to copy");
    break;
    
  case float_phase:
    if (DOSORT)
      fsort();
    for (i = 0; i < n; i++) {
      mpc_set_prec(mroot[i], DBL_MANT_DIG);
      mpc_set_cplx(mroot[i], froot[i]);
      rdpe_set_d(drad[i], frad[i]);
    }
    break;
    
  case dpe_phase:
    if (DOSORT)
      dsort();
    for (i = 0; i < n; i++) {
      mpc_set_prec(mroot[i], DBL_MANT_DIG);
      mpc_set_cdpe(mroot[i], droot[i]);
    }
    break;
    
  case mp_phase:
    if (DOSORT)
      msort();
    break;
    
  }
}

/*************************************************************
 *                     SUBROUTINE DUMP                       *
 *************************************************************/
void
dump(FILE * dmpstr)
{
  int i;
  
  fprintf(dmpstr, "\nDumping...\n");

  /* output current status */
  fprintf(dmpstr, 
	  "Phase=%d, In=%d, Out=%d, Uncertain=%d, Zero=%d, Clusters=%d\n",
	  lastphase, count[0], count[1], count[2], zero_roots, nclust);

  /* output current approximations */
  fprintf(dmpstr, "\nCurrent approximations:\n");  
  for (i = 0; i < n; i++) {
    fprintf(dmpstr, "%d:\t", i);

    switch (lastphase) {
    case no_phase:
    case float_phase:
      cplx_outln_str(dmpstr, froot[i]);
      break;
    
    case dpe_phase:
      cdpe_outln_str(dmpstr, droot[i]);
      break;
    
    case mp_phase:
      mpc_outln_str(dmpstr, 10, 0, mroot[i]);
      break;
    }
  }

  /* output radii */
  fprintf(dmpstr, "\nCurrent radii:\n");
  for (i = 0; i < n; i++) {
    fprintf(dmpstr, "%d:\t", i);

    switch (lastphase) {
    case no_phase:
    case float_phase:
      fprintf(dmpstr, "%e\n", frad[i]);
      break;
    
    case dpe_phase:
    case mp_phase:
      rdpe_outln_str(dmpstr, drad[i]);
      break;
    }
  }

  /* output position */
  fprintf(dmpstr, "\nPos:\t");
  for (i = 0; i < n; i++)
    fprintf(dmpstr, "%4d", i);

  /* output status information */
  fprintf(dmpstr, "\nStatus:\t");
  for (i = 0; i < n; i++)
    fprintf(dmpstr, "%4.3s", status[i]);

  /* output cluster information */
  fprintf(dmpstr, "\nClust:\t");
  for (i = 0; i < n; i++)
    fprintf(dmpstr, "%4d", clust[i]);

  fprintf(dmpstr, "\nPunt:\t");
  for (i = 0; i < nclust; i++)
    fprintf(dmpstr, "%4d", punt[i]);

  fprintf(dmpstr, "\n\n");
}

/*************************************************************
 *                     SUBROUTINE WARN                       *
 *************************************************************/
void
warn(char *s)
{
  if (DOWARN)
    fprintf(logstr, "%s", s);
}

/*************************************************************
 *                     SUBROUTINE VAERROR                    *
 *************************************************************/
void
error(int args, ...)
{
  va_list ap;
  char * s;
  
  warn("! ");		/* output error message */
  va_start(ap, args);
  while(args--) {
    s = va_arg(ap, char *);
    warn(s);		/* output error message */
  }
  va_end(ap);
  warn("\n");		/* output error message */

  dump(logstr);		/* dump status		*/
  exit(EXIT_FAILURE);	/* exit program         */
}
