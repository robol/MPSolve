/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */


#include <stdlib.h>
#include <math.h>
#include <mps/gmptools.h>

/**********************************************
*                  MPZ_T                      *
**********************************************/

#ifndef mpz_swap
void
mpz_swap (mpz_t z1, mpz_t z2)
/* z1 <-> z2 */
{
  mpz_t t;

  mpz_Move (t, z1);
  mpz_Move (z1, z2);
  mpz_Move (z2, t);
}
#endif

#ifndef mpz_tstbit
int
mpz_tstbit (mpz_t z, unsigned long int pos)
{
  size_t size;

  size = mpz_sizeinbase (z, 2);

  if (mpz_sgn (z) == 0 || size <= pos)
    return 0;

  if (mpz_scan1 (z, pos) == pos)
    return 1;
  else
    return 0;
}
#endif

/* vector support functions */

void
mpz_vinit (mpz_t v[], unsigned long int size)
{
  unsigned long int i;

  for (i = 0; i < size; i++)
    mpz_init (v[i]);
}

void
mpz_vclear (mpz_t v[], unsigned long int size)
{
  unsigned long int i;

  for (i = 0; i < size; i++)
    mpz_clear (v[i]);
}

/**********************************************
*                  MPQ_T                      *
**********************************************/

#ifndef mpq_swap
void
mpq_swap (mpq_t q1, mpq_t q2)
/* q1 <-> q2 */
{
  mpq_t t;

  mpq_Move (t, q1);
  mpq_Move (q1, q2);
  mpq_Move (q2, t);
}
#endif

#ifndef mpq_out_str
void
mpq_out_str (FILE * stream, int base, mpq_t q)
{
  mpz_out_str (stdout, base, mpq_numref (q));
  putc ('/', stream);
  mpz_out_str (stdout, base, mpq_denref (q));
}
#endif

/* vector support functions */

void
mpq_vinit (mpq_t v[], unsigned long int size)
{
  unsigned long int i;

  for (i = 0; i < size; i++)
    mpq_init (v[i]);
}

void
mpq_vclear (mpq_t v[], unsigned long int size)
{
  unsigned long int i;

  for (i = 0; i < size; i++)
    mpq_clear (v[i]);
}

/**********************************************
*                  MPF_T                      *
**********************************************/
#ifndef mpf_swap
void
mpf_swap (mpf_t f1, mpf_t f2)
/* f1 <-> f2 */
{
  mpf_t t;

  mpf_Move (t, f1);
  mpf_Move (f1, f2);
  mpf_Move (f2, t);
}
#endif

void
mpf_set_2dl (mpf_t f, double d, long int l)
{
  mpf_set_d (f, d);
  if (l >= 0)
    mpf_mul_2exp (f, f, l);
  else
    mpf_div_2exp (f, f, -l);
}

void
mpf_get_2dl (double *d, long int *l, mpf_t f)
{
  mp_exp_t e;
  double t;
  int i;

  /* pick mantissa and exponent from f */
  e = f->_mp_exp;
  f->_mp_exp = 0;
  t = mpf_get_d (f);
  f->_mp_exp = e;

  /* scale mantissa to (1/2, 1] */
  *d = frexp (t, &i);
  *l = e * mp_bits_per_limb + i;
}

long int
mpf_size_2 (mpf_t f)
{
  double d;
  long int l;

  /* pick mantissa and exponent from f */
  mpf_get_2dl (&d, &l, f);

  return l;
}

/* missing functions */

void
mpf_add_si (mpf_t r, mpf_t f, long int i)
{
  if (i >= 0)
    mpf_add_ui (r, f, i);
  else
    mpf_sub_ui (r, f, -i);
}

void
mpf_sub_si (mpf_t r, mpf_t f, long int i)
{
  if (i >= 0)
    mpf_sub_ui (r, f, i);
  else
    mpf_add_ui (r, f, -i);
}

void
mpf_si_sub (mpf_t r, long int i, mpf_t f)
{
  if (i >= 0)
    mpf_ui_sub (r, i, f);
  else
    {
      mpf_add_ui (r, f, -i);
      mpf_neg (r, r);
    }
}

void
mpf_mul_si (mpf_t r, mpf_t f, long int i)
{
  if (i >= 0)
    mpf_mul_ui (r, f, i);
  else
    {
      mpf_mul_ui (r, f, -i);
      mpf_neg (r, r);
    }
}

void
mpf_div_si (mpf_t r, mpf_t f, long int i)
{
  if (i >= 0)
    mpf_div_ui (r, f, i);
  else
    {
      mpf_div_ui (r, f, -i);
      mpf_neg (r, r);
    }
}

void
mpf_si_div (mpf_t r, long int i, mpf_t f)
{
  if (i >= 0)
    mpf_ui_div (r, i, f);
  else
    {
      mpf_ui_div (r, -i, f);
      mpf_neg (r, r);
    }
}

#ifndef mpf_pow_ui
void
mpf_pow_ui (mpf_t r, mpf_t f, register unsigned long int i)
/* r = f^i */
{
  mpf_t t;
  mpf_init2 (t, mpf_get_prec (r));

  mpf_set (t, f);

  if (i & 1)
    mpf_set (r, t);
  else
    mpf_set_ui (r, 1);
  i >>= 1;                      /* divide i by 2 */

  while (i)
    {
      mpf_sqr_eq (t);
      if (i & 1)
        mpf_mul_eq (r, t);
      i >>= 1;                  /* divide i by 2 */
    }

  mpf_clear (t);
}
#endif

void
mpf_pow_si (mpf_t r, mpf_t f, long int i)
/* r = f^i, i integer */
{
  if (i >= 0)
    mpf_pow_ui (r, f, i);
  else
    {
      mpf_pow_ui (r, f, -i);
      mpf_inv_eq (r);
    }
}

/* vector support functions */

void
mpf_vinit (mpf_t * v, unsigned long int size)
{
  unsigned long int i;

  for (i = 0; i < size; i++)
    mpf_init (v[i]);
}

void
mpf_vinit2 (mpf_t v[], unsigned long int size, unsigned long int prec)
{
  unsigned long int i;

  for (i = 0; i < size; i++)
    mpf_init2 (v[i], prec);
}

void
mpf_vclear (mpf_t v[], unsigned long int size)
{
  unsigned long int i;

  for (i = 0; i < size; i++)
    mpf_clear (v[i]);
}
