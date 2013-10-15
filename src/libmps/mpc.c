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
#include <mps/mpc.h>
#include <mps/gmptools.h>
#include <mps/mt.h>
#include <mps/tools.h>

/***********************************************************
**              functions for mpc_t                       **
***********************************************************/

/* constructors / destructors */
void
mpc_init (mpc_t c)
{
  printf ("mpc_init () called\n");
  abort ();
  mpf_init (mpc_Re (c));
  mpf_init (mpc_Im (c));
}

void
mpc_init2 (mpc_t c, unsigned long int prec)
{
  prec = (prec <= 2) ? 53 : prec;
  mpf_init2 (mpc_Re (c), prec);
  mpf_init2 (mpc_Im (c), prec);
}

void
mpc_clear (mpc_t c)
{
  mpf_clear (mpc_Re (c));
  mpf_clear (mpc_Im (c));
}

/* precision related functions */
void
mpc_set_prec (mpc_t c, unsigned long int prec)
{
  prec = (prec <= 2) ? 53 : prec;
  if (mpf_get_prec (mpc_Re (c)) < prec)
    {
      mpf_set_prec (mpc_Re (c), prec);
      mpf_set_prec (mpc_Im (c), prec);
    }
  else
    {
      mpf_set_prec (mpc_Re (c), prec);
      mpf_set_prec (mpc_Im (c), prec);
    }
}

unsigned long int
mpc_get_prec (mpc_t c)
{
  return mpf_get_prec (mpc_Re (c));
}

void
mpc_set_prec_raw (mpc_t c, unsigned long int prec)
{
  mpf_set_prec_raw (mpc_Re (c), prec);
  mpf_set_prec_raw (mpc_Im (c), prec);
}

/* assignment functions */
void
mpc_swap (mpc_t c1, mpc_t c2)
/* c1 <-> c2 */
{
  mpc_t t;

  mpc_Move (t, c1);
  mpc_Move (c1, c2);
  mpc_Move (c2, t);
}

void
mpc_set (mpc_t rc, mpc_t c)
{
  mpf_set (mpc_Re (rc), mpc_Re (c));
  mpf_set (mpc_Im (rc), mpc_Im (c));
}

void
mpc_set_ui (mpc_t c, unsigned long int ir, unsigned long int ii)
{
  mpf_set_ui (mpc_Re (c), ir);
  mpf_set_ui (mpc_Im (c), ii);
}

void
mpc_set_si (mpc_t c, signed long int ir, signed long int ii)
{
  mpf_set_si (mpc_Re (c), ir);
  mpf_set_si (mpc_Im (c), ii);
}

void
mpc_set_d (mpc_t c, double dr, double di)
{
  mpf_set_d (mpc_Re (c), dr);
  mpf_set_d (mpc_Im (c), di);
}

void
mpc_set_z (mpc_t c, mpz_t zr, mpz_t zi)
{
  mpf_set_z (mpc_Re (c), zr);
  mpf_set_z (mpc_Im (c), zi);
}

void
mpc_set_q (mpc_t c, mpq_t qr, mpq_t qi)
{
  mpf_set_q (mpc_Re (c), qr);
  mpf_set_q (mpc_Im (c), qi);
}

void
mpc_set_f (mpc_t c, mpf_t fr, mpf_t fi)
{
  mpf_set (mpc_Re (c), fr);
  mpf_set (mpc_Im (c), fi);
}

int
mpc_set_str (mpc_t c, char *sr, char *si, int base)
{
  if (!mpf_set_str (mpc_Re (c), sr, base))
    return -1;
  return mpf_set_str (mpc_Im (c), si, base);
}

void
mpc_init_set (mpc_t rc, mpc_t c)
{
  mpf_init_set (mpc_Re (rc), mpc_Re (c));
  mpf_init_set (mpc_Im (rc), mpc_Im (c));
}

void
mpc_init_set_ui (mpc_t c, unsigned long int ir, unsigned long int ii)
{
  mpf_init_set_ui (mpc_Re (c), ir);
  mpf_init_set_ui (mpc_Im (c), ii);
}

void
mpc_init_set_si (mpc_t c, signed long int ir, signed long int ii)
{
  mpf_init_set_si (mpc_Re (c), ir);
  mpf_init_set_si (mpc_Im (c), ii);
}

void
mpc_init_set_d (mpc_t c, double dr, double di)
{
  mpf_init_set_d (mpc_Re (c), dr);
  mpf_init_set_d (mpc_Im (c), di);
}

void
mpc_init_set_f (mpc_t c, mpf_t fr, mpf_t fi)
{
  mpf_init_set (mpc_Re (c), fr);
  mpf_init_set (mpc_Im (c), fi);
}

int
mpc_init_set_str (mpc_t c, char *sr, char *si, int base)
{
  if (!mpf_init_set_str (mpc_Re (c), sr, base))
    return -1;
  return mpf_init_set_str (mpc_Im (c), si, base);
}

/* unary functions */
void
mpc_neg (mpc_t rc, mpc_t c)
{
  mpf_neg (mpc_Re (rc), mpc_Re (c));
  mpf_neg (mpc_Im (rc), mpc_Im (c));
}

void
mpc_smod (mpf_t f, mpc_t c)
{
  mpf_t t;
  mpf_init2 (t, mpf_get_prec (f));

  mpf_mul (f, mpc_Re (c), mpc_Re (c));
  mpf_mul (t, mpc_Im (c), mpc_Im (c));
  mpf_add (f, f, t);

  mpf_clear (t);
}

void
mpc_rmod (rdpe_t r, mpc_t c)
{
  cdpe_t cdtmp;
  mpc_get_cdpe (cdtmp, c);
  cdpe_mod (r, cdtmp);
}

void
mpc_mod (mpf_t f, mpc_t c)
{
  mpf_t t;
  mpf_init2 (t, mpf_get_prec (f));

  mpf_mul (f, mpc_Re (c), mpc_Re (c));
  mpf_mul (t, mpc_Im (c), mpc_Im (c));
  mpf_add (f, f, t);
  mpf_sqrt (f, f);

  mpf_clear (t);
}

void
mpc_con (mpc_t rc, mpc_t c)
{
  mpc_set (rc, c);
  mpf_neg (mpc_Im (rc), mpc_Im (rc));
}

void
mpc_inv (mpc_t rc, mpc_t c)
{
  mpf_t f;
  mpf_init2 (f, mpf_get_prec (mpc_Re (rc)));

  mpc_smod (f, c);
  mpc_con (rc, c);
  mpf_div (mpc_Re (rc), mpc_Re (rc), f);
  mpf_div (mpc_Im (rc), mpc_Im (rc), f);

  mpf_clear (f);
}

void
mpc_inv2 (mpc_t rc, mpc_t c)
{
  mpf_t f;
  mpf_init2 (f, mpf_get_prec (mpc_Re (rc)));

  mpc_smod (f, c);
  mpf_ui_div (f, 1L, f);
  mpc_con (rc, c);
  mpf_mul (mpc_Re (rc), mpc_Re (rc), f);
  mpf_mul (mpc_Im (rc), mpc_Im (rc), f);

  mpf_clear (f);
}

void
mpc_sqr (mpc_t rc, mpc_t c)
{
  mpf_t f;
  mpf_init2 (f, mpf_get_prec (mpc_Re (rc)));

  mpf_mul (f, mpc_Re (c), mpc_Im (c));
  mpf_mul (mpc_Re (rc), mpc_Re (c), mpc_Re (c));
  mpf_mul (mpc_Im (rc), mpc_Im (c), mpc_Im (c));

  mpf_sub (mpc_Re (rc), mpc_Re (rc), mpc_Im (rc));
  mpf_mul_2exp (mpc_Im (rc), f, 1);

  mpf_clear (f);
}

void
mpc_rot (mpc_t rc, mpc_t c)
{
  mpf_t f;
  mpf_init2 (f, mpf_get_prec (mpc_Re (rc)));

  mpf_set (f, mpc_Re (c));
  mpf_set (mpc_Re (rc), mpc_Im (c));
  mpf_set (mpc_Im (rc), f);
  mpf_neg (mpc_Re (rc), mpc_Re (rc));

  mpf_clear (f);
}

void
mpc_flip (mpc_t rc, mpc_t c)
{
  mpf_t f;
  mpf_init2 (f, mpf_get_prec (mpc_Re (rc)));

  mpf_set (f, mpc_Re (c));
  mpf_set (mpc_Re (rc), mpc_Im (c));
  mpf_set (mpc_Im (rc), f);

  mpf_clear (f);
}

/* binary functions */
void
mpc_add (mpc_t rc, mpc_t c1, mpc_t c2)
{
  mpf_add (mpc_Re (rc), mpc_Re (c1), mpc_Re (c2));
  mpf_add (mpc_Im (rc), mpc_Im (c1), mpc_Im (c2));
}

void
mpc_add_f (mpc_t rc, mpc_t c, mpf_t f)
{
  mpf_add (mpc_Re (rc), mpc_Re (c), f);
}

void
mpc_add_ui (mpc_t rc, mpc_t c, unsigned long int r, unsigned long int i)
{
  mpf_add_ui (mpc_Re (rc), mpc_Re (c), r);
  mpf_add_ui (mpc_Im (rc), mpc_Im (c), i);
}

void
mpc_sub (mpc_t rc, mpc_t c1, mpc_t c2)
{
  mpf_sub (mpc_Re (rc), mpc_Re (c1), mpc_Re (c2));
  mpf_sub (mpc_Im (rc), mpc_Im (c1), mpc_Im (c2));
}

void
mpc_sub_f (mpc_t rc, mpc_t c, mpf_t f)
{
  mpf_sub (mpc_Re (rc), mpc_Re (c), f);
}

void
mpc_f_sub (mpc_t rc, mpf_t f, mpc_t c)
{
  mpf_sub (mpc_Re (rc), f, mpc_Re (c));
}

void
mpc_sub_ui (mpc_t rc, mpc_t c, unsigned long int r, unsigned long int i)
{
  mpf_sub_ui (mpc_Re (rc), mpc_Re (c), r);
  mpf_sub_ui (mpc_Im (rc), mpc_Im (c), i);
}

void
mpc_ui_sub (mpc_t rc, unsigned long int r, unsigned long int i, mpc_t c)
{
  mpf_ui_sub (mpc_Re (rc), r, mpc_Re (c));
  mpf_ui_sub (mpc_Im (rc), i, mpc_Im (c));
}

void
mpc_mul (mpc_t rc, mpc_t c1, mpc_t c2)
{
  mpf_t s1, s2, s3;
  unsigned long int i;

  i = mpf_get_prec (mpc_Re (rc));
  mpf_init2 (s1, i);
  mpf_init2 (s2, i);
  mpf_init2 (s3, i);

  mpf_set (s1, mpc_Re (c1));
  mpf_sub (s1, s1, mpc_Im (c1));
  mpf_set (s2, mpc_Re (c2));
  mpf_add (s2, s2, mpc_Im (c2));
  mpf_mul (s1, s1, s2);
  mpf_mul (s2, mpc_Re (c1), mpc_Im (c2));
  mpf_mul (s3, mpc_Im (c1), mpc_Re (c2));
  mpf_sub (mpc_Re (rc), s1, s2);
  mpf_add (mpc_Re (rc), mpc_Re (rc), s3);
  mpf_add (mpc_Im (rc), s2, s3);

  mpf_clear (s3);
  mpf_clear (s2);
  mpf_clear (s1);
}

void
mpc_mul_f (mpc_t rc, mpc_t c, mpf_t f)
{
  mpf_mul (mpc_Re (rc), mpc_Re (c), f);
  mpf_mul (mpc_Im (rc), mpc_Im (c), f);
}

void
mpc_mul_ui (mpc_t rc, mpc_t c, unsigned long int i)
{
  mpf_mul_ui (mpc_Re (rc), mpc_Re (c), i);
  mpf_mul_ui (mpc_Im (rc), mpc_Im (c), i);
}

void
mpc_mul_2exp (mpc_t rc, mpc_t c, unsigned long int i)
{
  mpf_mul_2exp (mpc_Re (rc), mpc_Re (c), i);
  mpf_mul_2exp (mpc_Im (rc), mpc_Im (c), i);
}

void
mpc_div (mpc_t rc, mpc_t c1, mpc_t c2)
{
  mpc_t t;
  mpc_init2 (t, mpf_get_prec (mpc_Re (rc)));

  mpc_inv (t, c2);
  mpc_mul (rc, c1, t);

  mpc_clear (t);
}

void
mpc_div_f (mpc_t rc, mpc_t c, mpf_t f)
{
  mpf_div (mpc_Re (rc), mpc_Re (c), f);
  mpf_div (mpc_Im (rc), mpc_Im (c), f);
}

void
mpc_f_div (mpc_t rc, mpf_t f, mpc_t c)
{
  mpc_t t;
  mpc_init2 (t, mpf_get_prec (mpc_Re (rc)));

  mpc_inv (t, c);
  mpc_mul_f (rc, t, f);

  mpc_clear (t);
}

void
mpc_div_ui (mpc_t rc, mpc_t c, unsigned long int i)
{
  mpf_div_ui (mpc_Re (rc), mpc_Re (c), i);
  mpf_div_ui (mpc_Im (rc), mpc_Im (c), i);
}

void
mpc_ui_div (mpc_t rc, unsigned long int i, mpc_t c)
{
  mpf_ui_div (mpc_Re (rc), i, mpc_Re (c));
  mpf_ui_div (mpc_Im (rc), i, mpc_Im (c));
}

void
mpc_div_2exp (mpc_t rc, mpc_t c, unsigned long int i)
{
  mpf_div_2exp (mpc_Re (rc), mpc_Re (c), i);
  mpf_div_2exp (mpc_Im (rc), mpc_Im (c), i);
}

void
mpc_pow_si (mpc_t rc, mpc_t c, register signed long int i)
/* rc = c^i, i integer */
{
  mpc_t t;
  mpc_init2 (t, mpf_get_prec (mpc_Re (rc)));

  mpc_set (t, c);
  if (i < 0)
    {
      mpc_inv_eq (t);
      i = -i;
    }
  if (i & 1)
    mpc_set (rc, t);
  else
    mpc_set_ui (rc, 1, 0);
  i >>= 1;                      /* divide i by 2 */

  while (i)
    {
      mpc_sqr_eq (t);
      if (i & 1)
        mpc_mul_eq (rc, t);
      i >>= 1;                  /* divide i by 2 */
    }

  mpc_clear (t);
}

/* op= style functions */
void
mpc_smod_eq (mpc_t c)
{
  mpf_mul (mpc_Re (c), mpc_Re (c), mpc_Re (c));
  mpf_mul (mpc_Im (c), mpc_Im (c), mpc_Im (c));
  mpf_add (mpc_Re (c), mpc_Re (c), mpc_Im (c));
  mpf_set_ui (mpc_Im (c), 0);
}

void
mpc_mod_eq (mpc_t c)
{
  mpf_mul (mpc_Re (c), mpc_Re (c), mpc_Re (c));
  mpf_mul (mpc_Im (c), mpc_Im (c), mpc_Im (c));
  mpf_add (mpc_Re (c), mpc_Re (c), mpc_Im (c));
  mpf_sqrt (mpc_Re (c), mpc_Re (c));
  mpf_set_ui (mpc_Im (c), 0);
}

void
mpc_rot_eq (mpc_t c)
{
  mpf_t f;

  mpf_Move (f, mpc_Re (c));
  mpf_Move (mpc_Re (c), mpc_Im (c));
  mpf_Move (mpc_Im (c), f);
  mpf_neg (mpc_Re (c), mpc_Re (c));
}

void
mpc_flip_eq (mpc_t c)
{
  mpf_t f;

  mpf_Move (f, mpc_Re (c));
  mpf_Move (mpc_Re (c), mpc_Im (c));
  mpf_Move (mpc_Im (c), f);
}

/* comparison functions */
int
mpc_eq (mpc_t c1, mpc_t c2, unsigned long int i)
{
  if (!mpf_eq (mpc_Re (c1), mpc_Re (c2), i))
    return -1;
  return mpf_eq (mpc_Im (c1), mpc_Im (c2), i);
}

int
mpc_eq_zero (mpc_t c)
{
  if (mpf_sgn (mpc_Re (c)) || mpf_sgn (mpc_Im (c)))
    return 0;
  else
    return 1;
}

int
mpc_eq_one (mpc_t c)
{
  if (mpf_sgn (mpc_Im (c)) || mpf_cmp_ui (mpc_Re (c), 1))
    return 0;
  else
    return 1;
}

/* I/O functions */
size_t
mpc_out_str_2u (FILE * f, int base, size_t n_digits_r,
                size_t n_digits_i, mpc_t c)
/* unformatted output */
{
  if (f == NULL)
    f = stdout;
  if (!mpf_out_str (f, base, n_digits_r, mpc_Re (c)))
    return 0;
  if (fprintf (f, " ") < 0)
    return 0;
  if (!mpf_out_str (f, base, n_digits_i, mpc_Im (c)))
    return 0;
  else
    return 1;
}

size_t
mpc_out_str_2 (FILE * f, int base, size_t n_digits_r,
               size_t n_digits_i, mpc_t c)
/* formatted output */
{
  if (f == NULL)
    f = stdout;
  if (fprintf (f, "(") < 0)
    return 0;
  if (!mpf_out_str (f, base, n_digits_r, mpc_Re (c)))
    return 0;
  if (fprintf (f, ", ") < 0)
    return 0;
  if (!mpf_out_str (f, base, n_digits_i, mpc_Im (c)))
    return 0;
  if (fprintf (f, ")") < 0)
    return 0;
  else
    return 1;
}

size_t
mpc_inp_str_u (mpc_t c, FILE * f, int base)
/* unformatted input */
{
  if (f == NULL)
    f = stdin;
  if (!mpf_inp_str (mpc_Re (c), f, base))
    return 0;
  if (fscanf (f, " ") < 0)
    return 0;
  if (!mpf_inp_str (mpc_Im (c), f, base))
    return 0;
  else
    return 1;
}

size_t
mpc_inp_str (mpc_t c, FILE * f, int base)
/* formatted input */
{
  if (f == NULL)
    f = stdin;
  if (fscanf (f, "(") < 0)
    return 0;
  if (!mpf_inp_str (mpc_Re (c), f, base))
    return 0;
  if (fscanf (f, ", ") < 0)
    return 0;
  if (!mpf_inp_str (mpc_Im (c), f, base))
    return 0;
  if (fscanf (f, ")") < 0)
    return 0;
  else
    return 1;
}

/* vector functions */
void
mpc_vinit (mpc_t v[], long size)
{
  long i;

  for (i = 0; i < size; i++)
    mpc_init (v[i]);
}

void
mpc_vinit2 (mpc_t v[], long size, long prec)
{
  long i;

  for (i = 0; i < size; i++)
    mpc_init2 (v[i], prec);
}

void
mpc_vclear (mpc_t v[], long size)
{
  long i;

  for (i = 0; i < size; i++)
    mpc_clear (v[i]);
  /* free(v); */
}

/***********************************************************
**                                                        **
***********************************************************/
