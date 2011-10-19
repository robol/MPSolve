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

#include <float.h>
#include <string.h>
#include <ctype.h>
#include <mps/gmptools.h>
#include <mps/poly.h>

/***********************************************************
*       SUBROUTINE READ_POLY                               *
***********************************************************/
void
mps_read_poly (mps_status * s, FILE * instr, mpspoly_t p)
{
  int indx, num_coeff, i, read_elements;

  /* skip blank or comment lines */
  while ((i = fgetc (instr)) == '!' || isspace (i))
    if (i == '!')
      while (fgetc (instr) != '\n')
        /* skip */ ;
  ungetc (i, instr);

  /* read data type */
  read_elements = fscanf (instr, "%3s", p->data_type);
  if (!read_elements)
    {
      mps_error (s, 1, "Error reading data type, aborting.");
    }

  /* read and convert prec_in from base 10 to base 2 */
  if (s->input_config->prec == -1)
    {
      /* Read input precision and abort on failure */
      read_elements = fscanf (instr, "%ld", &(p->prec_in));
      if (!read_elements)
        {
          mps_error (s, 1, "Error reading input precision, aborting.");
        }
      p->prec_in = (long) (p->prec_in * LOG2_10);
    }
  else
    {                           /* override input precision */
      read_elements = fscanf (instr, "%*d");
      p->prec_in = s->input_config->prec;  /* load default */
    }

  /* read degree */
  read_elements = fscanf (instr, "%d", &(p->deg));
  if (!read_elements)
    {
      mps_error (s, 1, "Error reading degree of polynomial, aborting.\n");
    }
  p->n = p->deg;

  /* if sparse read num of coeff */
  if (p->data_type[0] == 's')
    {
      read_elements = fscanf (instr, "%d", &num_coeff);
      if (!read_elements)
        {
          mps_error (s, 1, "Error reading number of coefficients, aborting.");
        }
    }
  else
    num_coeff = p->deg + 1;

  /* validate polynomial data */
  mps_validate_poly (s, p, num_coeff);

  /* setup floating point multiprecision */
  if (p->prec_in != 0)
    mps_mp_set_prec (s, p->prec_in);
  else
    mps_mp_set_prec (s, 2 * DBL_MANT_DIG);

  /* no need to read coefficients if user polynomial */
  if (p->data_type[0] == 'u')
    return;

  /* allocate polynomial vector */
  mps_allocate_poly (s, p);

  /* setup sparsity vector */
  if (p->data_type[0] == 's')
    for (i = 0; i <= p->deg + 1; i++)
      p->spar[i] = false;

  /* read coefficients */
  for (i = 0; i < num_coeff; i++)
    {

      if (p->data_type[0] == 's')
        {
          read_elements = fscanf (instr, "%d", &indx);
          if (!read_elements)
            {
              mps_error (s, 1,
                         "Error reading the polynomial coefficients, aborting.");
            }
          p->spar[indx] = true;
        }
      else
        indx = i;

      switch (p->data_type[1])
        {                       /* switch 1 */

        case 'r':              /* Real */

          switch (p->data_type[2])
            {                   /* switch 2 */

            case 'i':          /* Real - Integer Coefs */
              if (!mpz_inp_str (p->mip_r[indx], instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              break;

            case 'q':          /* Real - Rational Coefs */
              if (!mpz_inp_str (mpq_numref (p->mqp_r[indx]), instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              if (!mpz_inp_str (mpq_denref (p->mqp_r[indx]), instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              mpq_canonicalize (p->mqp_r[indx]);
              break;

            case 'f':          /* Real - Big/Float Coefs */
              if (!mpf_inp_str (p->mfpr[indx], instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              break;

            }                   /* switch 2 */
          break;

        case 'c':              /* Complex */

          switch (p->data_type[2])
            {                   /* switch 3 */

            case 'i':          /* Complex - Integer Coefs */
              if (!mpz_inp_str (p->mip_r[indx], instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              if (!mpz_inp_str (p->mip_i[indx], instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              break;

            case 'q':          /* Complex - Rational Coefs */
              if (!mpz_inp_str (mpq_numref (p->mqp_r[indx]), instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              if (!mpz_inp_str (mpq_denref (p->mqp_r[indx]), instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              mpq_canonicalize (p->mqp_r[indx]);
              if (!mpz_inp_str (mpq_numref (p->mqp_i[indx]), instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              if (!mpz_inp_str (mpq_denref (p->mqp_i[indx]), instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              mpq_canonicalize (p->mqp_i[indx]);
              break;

            case 'f':          /* Complex - Big/Float Coefs */
              if (!mpf_inp_str (mpc_Re (p->mfpc[indx]), instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              if (!mpf_inp_str (mpc_Im (p->mfpc[indx]), instr, 10))
                mps_error (s, 1, "Error while reading coefficient");
              break;

            }                   /* switch 3 */
        }                       /* switch 1 */
    }                           /* for */

  /* real begin */
  if (p->data_type[1] == 'r' && p->data_type[2] == 'f')
    for (indx = 0; indx <= p->deg; indx++)
      {
        if (p->data_type[0] == 's' && !p->spar[indx])
          continue;
        mpf_set (mpc_Re (p->mfpc[indx]), p->mfpr[indx]);
        mpf_set_ui (mpc_Im (p->mfpc[indx]), 0);
      }
  /* real end */
}

/*********************************************************
*                   PROCEDURE VALIDATE_POLY
**********************************************************/
void
mps_validate_poly (mps_status * s, mpspoly_t p, int num_coeff)
{
  /* check polynomial type */
  if (!p->data_type[0] || !strchr ("sdu", p->data_type[0]))
    mps_error (s, 1, "Illegal polynomial type");

  /* check ring type */
  if (!p->data_type[1] || !strchr ("rc", p->data_type[1]))
    mps_error (s, 1, "Illegal ring type");

  /* check coefficient type */
  if (!p->data_type[2] || !strchr ("iqfdb", p->data_type[2]))
    mps_error (s, 1, "Illegal coefficient type");

  /* check precision */
  if (p->prec_in < 0 || (p->prec_in == 0 && p->data_type[2] == 'f'))
    mps_error (s, 1, "Invalid input precision");

  /* check degree */
  if (p->n <= 0)
    mps_error (s, 1, "Invalid degree");

  /* check num coeff */
  if (num_coeff <= 0 || num_coeff > p->n + 1)
    mps_error (s, 1, "Invalid number of coefficients");

  /* log */
  if (s->DOLOG)
    {

      /* degree */
      fprintf (s->logstr, "Degree= %d\n", p->n);

      /* type */
      switch (p->data_type[0])
        {
        case 'u':
          fprintf (s->logstr, "User polynomial\n");
          break;
        case 's':
          fprintf (s->logstr, "Sparse polynomial\n");
          fprintf (s->logstr, "Num_coeff.= %d\n", num_coeff);
          break;
        case 'd':
          fprintf (s->logstr, "Dense polynomial\n");
          break;
        }

      /* ring */
      switch (p->data_type[1])
        {
        case 'r':
          fprintf (s->logstr, "Real case\n");
          break;
        case 'c':
          fprintf (s->logstr, "Complex polynomial\n");
          break;
        }

      /* coefficient */
      switch (p->data_type[2])
        {
        case 'i':
          fprintf (s->logstr, "Integer case\n");
          break;
        case 'q':
          fprintf (s->logstr, "Rational case\n");
          break;
        case 'f':
          fprintf (s->logstr, "Float case\n");
          break;
        case 'd':
          fprintf (s->logstr, "DPE case\n");
          break;
        case 'b':
          fprintf (s->logstr, "Bigfloat case\n");
          break;
        }
    }
}

/***********************************************************
*       SUBROUTINE SET_POLY                                *
***********************************************************/
void
mps_set_poly (mps_status * s, mpspoly_t p)
{
  s->n = p->n;
  s->deg = p->deg;
  s->data_type = p->data_type;
  s->input_config->prec = p->prec_in;
  s->spar = p->spar;
  s->fpr = p->fpr;
  s->fpc = p->fpc;
  s->dpr = p->dpr;
  s->dpc = p->dpc;
  s->mip_r = p->mip_r;
  s->mip_i = p->mip_i;
  s->mqp_r = p->mqp_r;
  s->mqp_i = p->mqp_i;
  s->mfpr = p->mfpr;
  s->mfpc = p->mfpc;
}

/***********************************************************
*       SUBROUTINE UPDATE_POLY                            *
***********************************************************/
void
mps_update_poly (mps_status * s, mpspoly_t p)
{
  p->prec_in = s->input_config->prec;
  p->n = s->n;
}

/********************************************************
 *      SUBROUTINE ALLOCATE_POLY                        *
 ********************************************************
  assumes that data_type, deg and prec_in are set
 *******************************************************/
/*DAFARE da perfezionare senza allocare cose inutili */
void
mps_allocate_poly (mps_status * s, mpspoly_t p)
{
  int i;
  if (s->DOLOG)
    fprintf (s->logstr, "Allocating polynomial\n");

  p->spar = mps_boolean_valloc (p->deg + 2);

  p->fpr = double_valloc (p->deg + 1);
  p->fpc = cplx_valloc (p->deg + 1);

  p->dpr = rdpe_valloc (p->deg + 1);
  p->dpc = cdpe_valloc (p->deg + 1);

  p->mip_r = mpz_valloc (p->deg + 1);
  p->mip_i = mpz_valloc (p->deg + 1);
  for (i = 0; i <= p->deg; i++)
    {
      mpz_init (p->mip_r[i]);
      mpz_init (p->mip_i[i]);
    }

  p->mqp_r = mpq_valloc (p->deg + 1);
  p->mqp_i = mpq_valloc (p->deg + 1);
  for (i = 0; i <= p->deg; i++)
    {
      mpq_init (p->mqp_r[i]);
      mpq_init (p->mqp_i[i]);
    }

  p->mfpr = mpf_valloc (p->deg + 1);
  for (i = 0; i <= p->deg; i++)
    mpf_init2 (p->mfpr[i], p->prec_in);
  p->mfpc = mpc_valloc (p->deg + 1);
  for (i = 0; i <= p->deg; i++)
    mpc_init2 (p->mfpc[i], p->prec_in);
}

/********************************************************
 *      SUBROUTINE FREE_POLY                            *
 *******************************************************/
void
mps_free_poly (mps_status * s, mpspoly_t p)
{
  int i;

  if (s->DOLOG)
    fprintf (s->logstr, "Unallocating polynomial\n");

  free (p->spar);

  free (p->fpr);
  cplx_vfree (p->fpc);

  rdpe_vfree (p->dpr);
  cdpe_vfree (p->dpc);

  for (i = 0; i <= p->deg; i++)
    {
      mpz_clear (p->mip_r[i]);
      mpz_clear (p->mip_i[i]);
    }
  free (p->mip_r);
  free (p->mip_i);

  for (i = 0; i <= p->deg; i++)
    {
      mpq_clear (p->mqp_r[i]);
      mpq_clear (p->mqp_i[i]);
    }
  free (p->mqp_r);
  free (p->mqp_i);

  for (i = 0; i <= p->deg; i++)
    mpf_clear (p->mfpr[i]);
  free (p->mfpr);

  for (i = 0; i <= p->deg; i++)
    mpc_clear (p->mfpc[i]);
  free (p->mfpc);
}

/********************************************************
 *      SUBROUTINE GET_ROOTS                            *
 *******************************************************/
void
mps_get_roots (mps_status * s, mpsroots_t r)
{
  r->lastphase = s->lastphase;
  r->count = s->count;
  r->zero_roots = s->zero_roots;
  r->status = s->status;
  r->froot = s->froot;
  r->droot = s->droot;
  r->mroot = s->mroot;
  r->frad = s->frad;
  r->drad = s->drad;
}
