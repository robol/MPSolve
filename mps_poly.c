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

#include <string.h>
#include <ctype.h>
#include "mps_poly.h"

/***********************************************************
*       SUBROUTINE READ_POLY                               *
***********************************************************/
void
read_poly(FILE *instr, mpspoly_t p)
{
  int indx, num_coeff, i;

  /* skip blank or comment lines */
  while ((i = fgetc(instr)) == '!' || isspace(i))
    if(i == '!')
      while (fgetc(instr) != '\n')
        /* skip */ ;
  ungetc(i, instr);

  /* read data type */
  fscanf(instr, "%3s", p->data_type);

  /* read and convert prec_in from base 10 to base 2 */
  if (prec_in == -1) {
    fscanf(instr, "%ld", &(p->prec_in));
    p->prec_in = (long) (p->prec_in * LOG2_10);
  } else {		/* override input precision */
    fscanf(instr, "%*d");
    p->prec_in = prec_in; /* load default */
  }

  /* read degree */
  fscanf(instr, "%d", &(p->deg));
  p->n = p->deg;

  /* if sparse read num of coeff */
  if (p->data_type[0] == 's')
    fscanf(instr, "%d", &num_coeff);
  else
    num_coeff = p->deg + 1;

  /* validate polynomial data */
  validate_poly(p, num_coeff);

  /* setup floating point multiprecision */
  if (p->prec_in != 0)
    mp_set_prec(p->prec_in);
  else
    mp_set_prec(2 * DBL_MANT_DIG);

  /* no need to read coefficients if user polynomial */
  if (p->data_type[0] == 'u')
    return;

  /* allocate polynomial vector */
  allocate_poly(p);

  /* setup sparsity vector */
  if (p->data_type[0] == 's')
    for (i = 0; i <= p->deg + 1; i++)
      p->spar[i] = false;
  
  /* read coefficients */
  for (i = 0; i < num_coeff; i++) {

    if (p->data_type[0] == 's') {
      fscanf(instr, "%d", &indx);
      p->spar[indx] = true;
    } else
      indx = i;

    switch (p->data_type[1]) {	/* switch 1 */

    case 'r':			/* Real */

      switch (p->data_type[2]) {	/* switch 2 */

      case 'i':		/* Real - Integer Coefs */
	if (!mpz_inp_str(p->mip_r[indx], instr, 10))
	  error(1, "Error while reading coefficient");
	break;

      case 'q':		/* Real - Rational Coefs */
	if (!mpz_inp_str(mpq_numref(p->mqp_r[indx]), instr, 10))
	  error(1, "Error while reading coefficient");
	if (!mpz_inp_str(mpq_denref(p->mqp_r[indx]), instr, 10))
	  error(1, "Error while reading coefficient");
	mpq_canonicalize(p->mqp_r[indx]);
	break;

      case 'f':		/* Real - Big/Float Coefs */
	if (!mpf_inp_str(p->mfpr[indx], instr, 10))
	  error(1, "Error while reading coefficient");
	break;

      }				/* switch 2 */
      break;

    case 'c':			/* Complex */

      switch (p->data_type[2]) {	/* switch 3 */

      case 'i':		/* Complex - Integer Coefs */
	if (!mpz_inp_str(p->mip_r[indx], instr, 10))
	  error(1, "Error while reading coefficient");
	if (!mpz_inp_str(p->mip_i[indx], instr, 10))
	  error(1, "Error while reading coefficient");
	break;

      case 'q':		/* Complex - Rational Coefs */
	if (!mpz_inp_str(mpq_numref(p->mqp_r[indx]), instr, 10))
	  error(1, "Error while reading coefficient");
	if (!mpz_inp_str(mpq_denref(p->mqp_r[indx]), instr, 10))
	  error(1, "Error while reading coefficient");
	mpq_canonicalize(p->mqp_r[indx]);
	if (!mpz_inp_str(mpq_numref(p->mqp_i[indx]), instr, 10))
	  error(1, "Error while reading coefficient");
	if (!mpz_inp_str(mpq_denref(p->mqp_i[indx]), instr, 10))
	  error(1, "Error while reading coefficient");
	mpq_canonicalize(p->mqp_i[indx]);
	break;

      case 'f':		/* Complex - Big/Float Coefs */
	if (!mpf_inp_str(mpc_Re(p->mfpc[indx]), instr, 10))
	  error(1, "Error while reading coefficient");
	if (!mpf_inp_str(mpc_Im(p->mfpc[indx]), instr, 10))
	  error(1, "Error while reading coefficient");
	break;

      }				/* switch 3 */
    }				/* switch 1 */
  }				/* for */

  /* real begin */
  if (p->data_type[1] == 'r' && p->data_type[2] == 'f')
    for (indx = 0; indx <= p->deg; indx++) {
      if (p->data_type[0] == 's' && !p->spar[indx])
	continue;
      mpf_set(mpc_Re(p->mfpc[indx]), p->mfpr[indx]);
      mpf_set_ui(mpc_Im(p->mfpc[indx]), 0);
    }
  /* real end */
}

/*********************************************************
*                   PROCEDURE VALIDATE_POLY
**********************************************************/
void
validate_poly(mpspoly_t p, int num_coeff)
{
  /* check polynomial type */
  if (!p->data_type[0] || !strchr("sdu", p->data_type[0]))
    error(1, "Illegal polynomial type");

  /* check ring type */
  if (!p->data_type[1] || !strchr("rc", p->data_type[1]))
    error(1, "Illegal ring type");

  /* check coefficient type */
  if (!p->data_type[2] || !strchr("iqfdb", p->data_type[2]))
    error(1, "Illegal coefficient type");

  /* check precision */
  if (p->prec_in < 0 || (p->prec_in == 0 && p->data_type[2] == 'f'))
    error(1, "Invalid input precision");

  /* check degree */
  if (p->n <= 0)
    error(1, "Invalid degree");

  /* check num coeff */
  if (num_coeff <= 0 || num_coeff > p->n + 1)
    error(1, "Invalid number of coefficients");

  /* log */
  if (DOLOG) {

    /* degree */
    fprintf(logstr, "Degree= %d\n", p->n);

    /* type */
    switch (p->data_type[0]) {
    case 'u':
      fprintf(logstr, "User polynomial\n");
      break;
    case 's':
      fprintf(logstr, "Sparse polynomial\n");
      fprintf(logstr, "Num_coeff.= %d\n", num_coeff);
      break;
    case 'd':
      fprintf(logstr, "Dense polynomial\n");
      break;
    }

    /* ring */
    switch (p->data_type[1]) {
    case 'r':
      fprintf(logstr, "Real case\n");
      break;
    case 'c':
      fprintf(logstr, "Complex polynomial\n");
      break;
    }

    /* coefficient */
    switch (p->data_type[2]) {
    case 'i':
      fprintf(logstr, "Integer case\n");
      break;
    case 'q':
      fprintf(logstr, "Rational case\n");
      break;
    case 'f':
      fprintf(logstr, "Float case\n");
      break;
    case 'd':
      fprintf(logstr, "DPE case\n");
      break;
    case 'b':
      fprintf(logstr, "Bigfloat case\n");
      break;
    }
  }
}

/***********************************************************
*       SUBROUTINE SET_POLY                                *
***********************************************************/
void
set_poly(mpspoly_t p)
{
  n = p->n;  
  deg = p->deg;
  data_type = p->data_type;
  prec_in = p->prec_in;
  spar = p->spar;
  fpr = p->fpr;
  fpc = p->fpc;
  dpr = p->dpr;
  dpc = p->dpc;
  mip_r = p->mip_r;
  mip_i = p->mip_i;
  mqp_r = p->mqp_r;
  mqp_i = p->mqp_i;
  mfpr = p->mfpr;
  mfpc = p->mfpc;	
}

/***********************************************************
*       SUBROUTINE UPDATE_POLY                            *
***********************************************************/
void
update_poly(mpspoly_t p)
{
  p->prec_in = prec_in;
  p->n = n;  
}

/********************************************************
 *      SUBROUTINE ALLOCATE_POLY                        *
 ********************************************************
  assumes that data_type, deg and prec_in are set
 *******************************************************/
/*DAFARE da perfezionare senza allocare cose inutili */
void
allocate_poly(mpspoly_t p)
{
  int i;
  if (DOLOG)
    fprintf(logstr, "Allocating polynomial\n");

  p->spar = boolean_valloc(p->deg + 2);

  p->fpr = double_valloc(p->deg + 1);
  p->fpc = cplx_valloc(p->deg + 1);

  p->dpr = rdpe_valloc(p->deg + 1);
  p->dpc = cdpe_valloc(p->deg + 1);

  p->mip_r = mpz_valloc(p->deg + 1);
  p->mip_i = mpz_valloc(p->deg + 1);
  for (i = 0; i <= p->deg; i++) {
    mpz_init(p->mip_r[i]);
    mpz_init(p->mip_i[i]);
  }

  p->mqp_r = mpq_valloc(p->deg + 1);
  p->mqp_i = mpq_valloc(p->deg + 1);
  for (i = 0; i <= p->deg; i++) {
    mpq_init(p->mqp_r[i]);
    mpq_init(p->mqp_i[i]);
  }

  p->mfpr = mpf_valloc(p->deg + 1);
  for (i = 0; i <= p->deg; i++)
    mpf_init2(p->mfpr[i], p->prec_in);
  p->mfpc = mpc_valloc(p->deg + 1);
  for (i = 0; i <= p->deg; i++)
    mpc_init2(p->mfpc[i], p->prec_in);
}

/********************************************************
 *      SUBROUTINE FREE_POLY                            *
 *******************************************************/
void
free_poly(mpspoly_t p)
{
  int i;

  if (DOLOG)
    fprintf(logstr, "Unallocating polynomial\n");

  free(p->spar);

  free(p->fpr);
  cplx_vfree(p->fpc);

  rdpe_vfree(p->dpr);
  cdpe_vfree(p->dpc);

  for (i = 0; i <= p->deg; i++) {
    mpz_clear(p->mip_r[i]);
    mpz_clear(p->mip_i[i]);
  }
  free(p->mip_r);
  free(p->mip_i);

  for (i = 0; i <= p->deg; i++) {
    mpq_clear(p->mqp_r[i]);
    mpq_clear(p->mqp_i[i]);
  }
  free(p->mqp_r);
  free(p->mqp_i);

  for (i = 0; i <= p->deg; i++)
    mpf_clear(p->mfpr[i]);
  free(p->mfpr);

  for (i = 0; i <= p->deg; i++)
    mpc_clear(p->mfpc[i]);
  free(p->mfpc);
}

/********************************************************
 *      SUBROUTINE GET_ROOTS                            *
 *******************************************************/
void get_roots(mpsroots_t r)
{
  r->lastphase = lastphase;
  r->count = count;
  r->zero_roots = zero_roots;
  r->status = status;
  r->froot = froot;
  r->droot = droot;
  r->mroot = mroot;
  r->frad = frad;
  r->drad = drad;
}
