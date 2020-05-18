/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#include <stdlib.h>
#include <mps/mps.h>

#define MPS_MPF_TEMP_SIZE 6

struct mps_tls {
  pthread_t thread;
  mpf_t *data;
  long int precision;
  struct mps_tls *next;
};

typedef struct mps_tls mps_tls;

static pthread_once_t once_key_created = PTHREAD_ONCE_INIT;
static pthread_key_t key;

static void
mps_mpcf_cache_cleanup (void * pointer)
{
  mps_tls *ptr = pointer;
  int i;

  for (i = 0; i < MPS_MPF_TEMP_SIZE; i++)
    mpf_clear (ptr->data[i]);

  free(ptr->data);
  free (ptr);
}

static mps_tls *
create_new_mps_tls (long int precision_needed)
{
  mps_tls *ptr;
  int i;

  ptr = mps_new (mps_tls);

  ptr->data = mps_newv (mpf_t, MPS_MPF_TEMP_SIZE);
  ptr->precision = precision_needed;

  for (i = 0; i < MPS_MPF_TEMP_SIZE; i++)
    mpf_init2 (ptr->data[i], precision_needed);

  /* Set up a destructor for this data in case the thread exits */
  pthread_setspecific (key, ptr);

  return ptr;
}

static void
adjust_mps_tls_precision (mps_tls *ptr, long int precision_needed)
{
  int i;

  for (i = 0; i < MPS_MPF_TEMP_SIZE; i++)
    mpf_set_prec (ptr->data[i], precision_needed);

  ptr->precision = precision_needed;
}

static void
create_key (void)
{
  pthread_key_create (&key, mps_mpcf_cache_cleanup);
}

static mpf_t*
init (long int precision_needed)
{
  pthread_once (&once_key_created, create_key);
  mps_tls *ptr = pthread_getspecific (key);

  /* This means that we have to create a new entry */
  if (ptr == NULL)
    {
      ptr = create_new_mps_tls (precision_needed);
    }
  else if (ptr->precision < precision_needed || precision_needed < .25 * ptr->precision)
    {
      adjust_mps_tls_precision (ptr, precision_needed);
    }

  return ptr->data;
}


/***********************************************************
**              functions for mpcf_t                       **
***********************************************************/

/* constructors / destructors */
void
mpcf_init (mpcf_t c)
{
  printf ("mpcf_init () called\n");
  abort ();
  mpf_init (mpcf_Re (c));
  mpf_init (mpcf_Im (c));
}

void
mpcf_init2 (mpcf_t c, unsigned long int prec)
{
  prec = (prec <= 2) ? 53 : prec;
  mpf_init2 (mpcf_Re (c), prec);
  mpf_init2 (mpcf_Im (c), prec);
}

void
mpcf_clear (mpcf_t c)
{
  mpf_clear (mpcf_Re (c));
  mpf_clear (mpcf_Im (c));
}

/* precision related functions */
void
mpcf_set_prec (mpcf_t c, unsigned long int prec)
{
  prec = (prec <= 2) ? 53 : prec;
  if (mpf_get_prec (mpcf_Re (c)) < prec)
    {
      mpf_set_prec (mpcf_Re (c), prec);
      mpf_set_prec (mpcf_Im (c), prec);
    }
  else
    {
      mpf_set_prec (mpcf_Re (c), prec);
      mpf_set_prec (mpcf_Im (c), prec);
    }
}

unsigned long int
mpcf_get_prec (const mpcf_t c)
{
  return mpf_get_prec (mpcf_Re (c));
}

void
mpcf_set_prec_raw (mpcf_t c, unsigned long int prec)
{
  mpf_set_prec_raw (mpcf_Re (c), prec);
  mpf_set_prec_raw (mpcf_Im (c), prec);
}

/* assignment functions */
void
mpcf_swap (mpcf_t c1, mpcf_t c2)
/* c1 <-> c2 */
{
  mpcf_t t;

  mpcf_Move (t, c1);
  mpcf_Move (c1, c2);
  mpcf_Move (c2, t);
}

void
mpcf_set (mpcf_t rc, const mpcf_t c)
{
  mpf_set (mpcf_Re (rc), mpcf_Re (c));
  mpf_set (mpcf_Im (rc), mpcf_Im (c));
}

void
mpcf_set_ui (mpcf_t c, unsigned long int ir, unsigned long int ii)
{
  mpf_set_ui (mpcf_Re (c), ir);
  mpf_set_ui (mpcf_Im (c), ii);
}

void
mpcf_set_si (mpcf_t c, signed long int ir, signed long int ii)
{
  mpf_set_si (mpcf_Re (c), ir);
  mpf_set_si (mpcf_Im (c), ii);
}

void
mpcf_set_d (mpcf_t c, double dr, double di)
{
  mpf_set_d (mpcf_Re (c), dr);
  mpf_set_d (mpcf_Im (c), di);
}

void
mpcf_set_z (mpcf_t c, mpz_t zr, mpz_t zi)
{
  mpf_set_z (mpcf_Re (c), zr);
  mpf_set_z (mpcf_Im (c), zi);
}

void
mpcf_set_q (mpcf_t c, mpq_t qr, mpq_t qi)
{
  mpf_set_q (mpcf_Re (c), qr);
  mpf_set_q (mpcf_Im (c), qi);
}

void
mpcf_set_f (mpcf_t c, mpf_t fr, mpf_t fi)
{
  mpf_set (mpcf_Re (c), fr);
  mpf_set (mpcf_Im (c), fi);
}

int
mpcf_set_str (mpcf_t c, char *sr, char *si, int base)
{
  if (!mpf_set_str (mpcf_Re (c), sr, base))
    return -1;
  return mpf_set_str (mpcf_Im (c), si, base);
}

void
mpcf_init_set (mpcf_t rc, mpcf_t c)
{
  mpf_init_set (mpcf_Re (rc), mpcf_Re (c));
  mpf_init_set (mpcf_Im (rc), mpcf_Im (c));
}

void
mpcf_init_set_ui (mpcf_t c, unsigned long int ir, unsigned long int ii)
{
  mpf_init_set_ui (mpcf_Re (c), ir);
  mpf_init_set_ui (mpcf_Im (c), ii);
}

void
mpcf_init_set_si (mpcf_t c, signed long int ir, signed long int ii)
{
  mpf_init_set_si (mpcf_Re (c), ir);
  mpf_init_set_si (mpcf_Im (c), ii);
}

void
mpcf_init_set_d (mpcf_t c, double dr, double di)
{
  mpf_init_set_d (mpcf_Re (c), dr);
  mpf_init_set_d (mpcf_Im (c), di);
}

void
mpcf_init_set_f (mpcf_t c, mpf_t fr, mpf_t fi)
{
  mpf_init_set (mpcf_Re (c), fr);
  mpf_init_set (mpcf_Im (c), fi);
}

int
mpcf_init_set_str (mpcf_t c, char *sr, char *si, int base)
{
  if (!mpf_init_set_str (mpcf_Re (c), sr, base))
    return -1;
  return mpf_init_set_str (mpcf_Im (c), si, base);
}

/* unary functions */
void
mpcf_neg (mpcf_t rc, mpcf_t c)
{
  mpf_neg (mpcf_Re (rc), mpcf_Re (c));
  mpf_neg (mpcf_Im (rc), mpcf_Im (c));
}

void
mpcf_smod (mpf_t f, mpcf_t c)
{
  mpf_t *t = init (mpf_get_prec (f));

  mpf_mul (f, mpcf_Re (c), mpcf_Re (c));
  mpf_mul (*t, mpcf_Im (c), mpcf_Im (c));
  mpf_add (f, f, *t);
}

void
mpcf_rmod (rdpe_t r, mpcf_t c)
{
  cdpe_t cdtmp;

  mpcf_get_cdpe (cdtmp, c);
  cdpe_mod (r, cdtmp);
}

void
mpcf_mod (mpf_t f, mpcf_t c)
{
  mpf_t *t = init (mpf_get_prec (f));

  mpf_mul (f, mpcf_Re (c), mpcf_Re (c));
  mpf_mul (*t, mpcf_Im (c), mpcf_Im (c));
  mpf_add (f, f, *t);
  mpf_sqrt (f, f);
}

void
mpcf_con (mpcf_t rc, mpcf_t c)
{
  mpcf_set (rc, c);
  mpf_neg (mpcf_Im (rc), mpcf_Im (rc));
}

void
mpcf_inv (mpcf_t rc, mpcf_t c)
{
  mpf_t *f = init (mpf_get_prec (mpcf_Re (rc))) + 5;

  mpcf_smod (*f, c);
  mpcf_con (rc, c);
  mpf_div (mpcf_Re (rc), mpcf_Re (rc), *f);
  mpf_div (mpcf_Im (rc), mpcf_Im (rc), *f);
}

void
mpcf_inv2 (mpcf_t rc, mpcf_t c)
{
  mpf_t *f = init (mpf_get_prec (mpcf_Re (rc))) + 5;

  mpcf_smod (*f, c);
  mpf_ui_div (*f, 1L, *f);
  mpcf_con (rc, c);
  mpf_mul (mpcf_Re (rc), mpcf_Re (rc), *f);
  mpf_mul (mpcf_Im (rc), mpcf_Im (rc), *f);
}

void
mpcf_sqr (mpcf_t rc, mpcf_t c)
{
  mpf_t *f = init (mpf_get_prec (mpcf_Re (rc)));

  mpf_mul (*f, mpcf_Re (c), mpcf_Im (c));
  mpf_mul (mpcf_Re (rc), mpcf_Re (c), mpcf_Re (c));
  mpf_mul (mpcf_Im (rc), mpcf_Im (c), mpcf_Im (c));

  mpf_sub (mpcf_Re (rc), mpcf_Re (rc), mpcf_Im (rc));
  mpf_mul_2exp (mpcf_Im (rc), *f, 1);
}

void
mpcf_rot (mpcf_t rc, mpcf_t c)
{
  mpf_t *f = init (mpf_get_prec (mpcf_Re (rc)));

  mpf_set (*f, mpcf_Re (c));
  mpf_set (mpcf_Re (rc), mpcf_Im (c));
  mpf_set (mpcf_Im (rc), *f);
  mpf_neg (mpcf_Re (rc), mpcf_Re (rc));
}

void
mpcf_flip (mpcf_t rc, mpcf_t c)
{
  mpf_t f;

  mpf_init2 (f, mpf_get_prec (mpcf_Re (rc)));

  mpf_set (f, mpcf_Re (c));
  mpf_set (mpcf_Re (rc), mpcf_Im (c));
  mpf_set (mpcf_Im (rc), f);

  mpf_clear (f);
}

/* binary functions */
void
mpcf_add (mpcf_t rc, mpcf_t c1, mpcf_t c2)
{
  mpf_add (mpcf_Re (rc), mpcf_Re (c1), mpcf_Re (c2));
  mpf_add (mpcf_Im (rc), mpcf_Im (c1), mpcf_Im (c2));
}

void
mpcf_add_f (mpcf_t rc, mpcf_t c, mpf_t f)
{
  mpf_add (mpcf_Re (rc), mpcf_Re (c), f);
}

void
mpcf_add_ui (mpcf_t rc, mpcf_t c, unsigned long int r, unsigned long int i)
{
  mpf_add_ui (mpcf_Re (rc), mpcf_Re (c), r);
  mpf_add_ui (mpcf_Im (rc), mpcf_Im (c), i);
}

void
mpcf_sub (mpcf_t rc, mpcf_t c1, mpcf_t c2)
{
  mpf_sub (mpcf_Re (rc), mpcf_Re (c1), mpcf_Re (c2));
  mpf_sub (mpcf_Im (rc), mpcf_Im (c1), mpcf_Im (c2));
}

void
mpcf_sub_f (mpcf_t rc, mpcf_t c, mpf_t f)
{
  mpf_sub (mpcf_Re (rc), mpcf_Re (c), f);
}

void
mpcf_f_sub (mpcf_t rc, mpf_t f, mpcf_t c)
{
  mpf_sub (mpcf_Re (rc), f, mpcf_Re (c));
}

void
mpcf_sub_ui (mpcf_t rc, mpcf_t c, unsigned long int r, unsigned long int i)
{
  mpf_sub_ui (mpcf_Re (rc), mpcf_Re (c), r);
  mpf_sub_ui (mpcf_Im (rc), mpcf_Im (c), i);
}

void
mpcf_ui_sub (mpcf_t rc, unsigned long int r, unsigned long int i, mpcf_t c)
{
  mpf_ui_sub (mpcf_Re (rc), r, mpcf_Re (c));
  mpf_ui_sub (mpcf_Im (rc), i, mpcf_Im (c));
}

void
mpcf_mul (mpcf_t rc, mpcf_t c1, mpcf_t c2)
{
  mpf_t *s1, *s2, *s3;

  s1 = init (mpf_get_prec (mpcf_Re (rc)));
  s2 = s1 + 1;
  s3 = s2 + 1;

  mpf_sub (*s1, mpcf_Re (c1), mpcf_Im (c1));
  mpf_add (*s2, mpcf_Re (c2), mpcf_Im (c2));
  mpf_mul (*s1, *s1, *s2);
  mpf_mul (*s2, mpcf_Re (c1), mpcf_Im (c2));
  mpf_mul (*s3, mpcf_Im (c1), mpcf_Re (c2));
  mpf_sub (mpcf_Re (rc), *s1, *s2);
  mpf_add (mpcf_Re (rc), mpcf_Re (rc), *s3);
  mpf_add (mpcf_Im (rc), *s2, *s3);
}

void
mpcf_mul_f (mpcf_t rc, mpcf_t c, mpf_t f)
{
  mpf_mul (mpcf_Re (rc), mpcf_Re (c), f);
  mpf_mul (mpcf_Im (rc), mpcf_Im (c), f);
}

void
mpcf_mul_ui (mpcf_t rc, mpcf_t c, unsigned long int i)
{
  mpf_mul_ui (mpcf_Re (rc), mpcf_Re (c), i);
  mpf_mul_ui (mpcf_Im (rc), mpcf_Im (c), i);
}

void
mpcf_mul_2exp (mpcf_t rc, mpcf_t c, unsigned long int i)
{
  mpf_mul_2exp (mpcf_Re (rc), mpcf_Re (c), i);
  mpf_mul_2exp (mpcf_Im (rc), mpcf_Im (c), i);
}

void
mpcf_div (mpcf_t rc, mpcf_t c1, mpcf_t c2)
{
  /* Find a good offset so that we don't pollute mul and inv registers */
  mpcf_t t;

  mpcf_init2 (t, mpf_get_prec (mpcf_Re (rc)));

  mpcf_inv (t, c2);
  mpcf_mul (rc, c1, t);
  mpcf_clear (t);
}

void
mpcf_div_f (mpcf_t rc, mpcf_t c, mpf_t f)
{
  mpf_div (mpcf_Re (rc), mpcf_Re (c), f);
  mpf_div (mpcf_Im (rc), mpcf_Im (c), f);
}

void
mpcf_f_div (mpcf_t rc, mpf_t f, mpcf_t c)
{
  mpcf_t t;

  mpcf_init2 (t, mpf_get_prec (mpcf_Re (rc)));

  mpcf_inv (t, c);
  mpcf_mul_f (rc, t, f);

  mpcf_clear (t);
}

void
mpcf_div_ui (mpcf_t rc, mpcf_t c, unsigned long int i)
{
  mpf_div_ui (mpcf_Re (rc), mpcf_Re (c), i);
  mpf_div_ui (mpcf_Im (rc), mpcf_Im (c), i);
}

void
mpcf_ui_div (mpcf_t rc, unsigned long int i, mpcf_t c)
{
  mpcf_set_ui(rc, i, 0L);
  mpcf_div_eq(rc, c);
}

void
mpcf_div_2exp (mpcf_t rc, mpcf_t c, unsigned long int i)
{
  mpf_div_2exp (mpcf_Re (rc), mpcf_Re (c), i);
  mpf_div_2exp (mpcf_Im (rc), mpcf_Im (c), i);
}

void
mpcf_pow_si (mpcf_t rc, mpcf_t c, register signed long int i)
/* rc = c^i, i integer */
{
  mpcf_t t;

  mpcf_init2 (t, mpf_get_prec (mpcf_Re (rc)));

  mpcf_set (t, c);
  if (i < 0)
    {
      mpcf_inv_eq (t);
      i = -i;
    }
  if (i & 1)
    mpcf_set (rc, t);
  else
    mpcf_set_ui (rc, 1, 0);
  i >>= 1;                      /* divide i by 2 */

  while (i)
    {
      mpcf_sqr_eq (t);
      if (i & 1)
        mpcf_mul_eq (rc, t);
      i >>= 1;                  /* divide i by 2 */
    }

  mpcf_clear (t);
}

/* op= style functions */
void
mpcf_smod_eq (mpcf_t c)
{
  mpf_mul (mpcf_Re (c), mpcf_Re (c), mpcf_Re (c));
  mpf_mul (mpcf_Im (c), mpcf_Im (c), mpcf_Im (c));
  mpf_add (mpcf_Re (c), mpcf_Re (c), mpcf_Im (c));
  mpf_set_ui (mpcf_Im (c), 0);
}

void
mpcf_mod_eq (mpcf_t c)
{
  mpf_mul (mpcf_Re (c), mpcf_Re (c), mpcf_Re (c));
  mpf_mul (mpcf_Im (c), mpcf_Im (c), mpcf_Im (c));
  mpf_add (mpcf_Re (c), mpcf_Re (c), mpcf_Im (c));
  mpf_sqrt (mpcf_Re (c), mpcf_Re (c));
  mpf_set_ui (mpcf_Im (c), 0);
}

void
mpcf_rot_eq (mpcf_t c)
{
  mpf_t f;

  mpf_Move (f, mpcf_Re (c));
  mpf_Move (mpcf_Re (c), mpcf_Im (c));
  mpf_Move (mpcf_Im (c), f);
  mpf_neg (mpcf_Re (c), mpcf_Re (c));
}

void
mpcf_flip_eq (mpcf_t c)
{
  mpf_t f;

  mpf_Move (f, mpcf_Re (c));
  mpf_Move (mpcf_Re (c), mpcf_Im (c));
  mpf_Move (mpcf_Im (c), f);
}

/* comparison functions */
int
mpcf_eq (mpcf_t c1, mpcf_t c2, unsigned long int i)
{
  if (!mpf_eq (mpcf_Re (c1), mpcf_Re (c2), i))
    return -1;
  return mpf_eq (mpcf_Im (c1), mpcf_Im (c2), i);
}

int
mpcf_eq_zero (mpcf_t c)
{
  if (mpf_sgn (mpcf_Re (c)) || mpf_sgn (mpcf_Im (c)))
    return 0;
  else
    return 1;
}

int
mpcf_eq_one (mpcf_t c)
{
  if (mpf_sgn (mpcf_Im (c)) || mpf_cmp_ui (mpcf_Re (c), 1))
    return 0;
  else
    return 1;
}

/* I/O functions */
size_t
mpcf_out_str_2u (FILE * f, int base, size_t n_digits_r,
                size_t n_digits_i, mpcf_t c)
/* unformatted output */
{
  if (f == NULL)
    f = stdout;
  if (!mpf_out_str (f, base, n_digits_r, mpcf_Re (c)))
    return 0;
  if (fprintf (f, " ") < 0)
    return 0;
  if (!mpf_out_str (f, base, n_digits_i, mpcf_Im (c)))
    return 0;
  else
    return 1;
}

size_t
mpcf_out_str_2 (FILE * f, int base, size_t n_digits_r,
               size_t n_digits_i, mpcf_t c)
/* formatted output */
{
  if (f == NULL)
    f = stdout;
  if (fprintf (f, "(") < 0)
    return 0;
  if (!mpf_out_str (f, base, n_digits_r, mpcf_Re (c)))
    return 0;
  if (fprintf (f, ", ") < 0)
    return 0;
  if (!mpf_out_str (f, base, n_digits_i, mpcf_Im (c)))
    return 0;
  if (fprintf (f, ")") < 0)
    return 0;
  else
    return 1;
}

size_t
mpcf_inp_str_u (mpcf_t c, FILE * f, int base)
/* unformatted input */
{
  if (f == NULL)
    f = stdin;
  if (!mpf_inp_str (mpcf_Re (c), f, base))
    return 0;
  if (fscanf (f, " ") < 0)
    return 0;
  if (!mpf_inp_str (mpcf_Im (c), f, base))
    return 0;
  else
    return 1;
}

size_t
mpcf_inp_str (mpcf_t c, FILE * f, int base)
/* formatted input */
{
  if (f == NULL)
    f = stdin;
  if (fscanf (f, "(") < 0)
    return 0;
  if (!mpf_inp_str (mpcf_Re (c), f, base))
    return 0; 
  if (fscanf (f, ", ") < 0)
    return 0;
  if (!mpf_inp_str (mpcf_Im (c), f, base)) 
    return 0; 
  if (fscanf (f, ")") < 0)
    return 0;
  else
    return 1;
}

/* vector functions */
void
mpcf_vinit (mpcf_t v[], long size)
{
  long i;

  for (i = 0; i < size; i++)
    mpcf_init (v[i]);
}

void
mpcf_vinit2 (mpcf_t v[], long size, long prec)
{
  long i;

  for (i = 0; i < size; i++)
    mpcf_init2 (v[i], prec);
}

void
mpcf_vclear (mpcf_t v[], long size)
{
  long i;

  for (i = 0; i < size; i++)
    mpcf_clear (v[i]);
  /* free(v); */
}

/***********************************************************
**                                                        **
***********************************************************/
