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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <mps/mt.h>
#include <mps/mps.h>

/***********************************************************
**              machine dependent constants               **
***********************************************************/

#define NBT           DBL_MANT_DIG

/***********************************************************
**              input/output formats                      **
***********************************************************/

#define FMTE           "% 18.14e"
#define FMTF           "% 16.14f"
#define CPLX_OUT_UFMT  FMTE " " FMTE
#define CPLX_OUT_FMT   "(" FMTE ", " FMTE ")"
#define CPLX_INP_UFMT  "%le %le"
#define CPLX_INP_FMT   "(%le, %le)"
#define RDPE_OUT_UFMT  FMTF "e%+04li"
#define RDPE_OUT_FMT   FMTF "x%+04li"
#define RDPE_INP_UFMT  "%lf e %ld"
#define RDPE_INP_FMT   "%lf x %ld"
#define CDPE_OUT_UFMT  RDPE_OUT_UFMT " " RDPE_OUT_UFMT
#define CDPE_OUT_FMT   "(" RDPE_OUT_FMT ", " RDPE_OUT_FMT ")"
#define CDPE_INP_UFMT  RDPE_INP_UFMT " " RDPE_INP_UFMT
#define CDPE_INP_FMT   "(" RDPE_INP_FMT ", " RDPE_INP_FMT ")"
#define DEF_STR_SIZE   80

/***********************************************************
**              numerical constants                       **
***********************************************************/

#define LOG10_2      0.30102999566398119521
#define LOG2_10      3.32192809488736234787
#define LOG_2        0.69314718055994530941
#define LOG_10       2.30258509299404568401

/***********************************************************
**              functions for cplx_t                      **
***********************************************************/

#ifdef MPS_USE_BUILTIN_COMPLEX

/* base constants */
const cplx_t cplx_zero = { {0.0, 0.0} };
const cplx_t cplx_one = { {1.0, 0.0} };
const cplx_t cplx_i = { {0.0, 1.0} };

void
cplx_d (cplx_t temp_cplx, double r, double i)
/* return (r, i) */
{
  cplx_Re (temp_cplx) = r;
  cplx_Im (temp_cplx) = i;
}

void
cplx_clear (cplx_t x)
/* x = 0 + I 0 */
{
  cplx_Move (x, cplx_zero);
}

void
cplx_set (cplx_t rx, const cplx_t x)
/* rx = x */
{
  cplx_Move (rx, x);
}

void
cplx_set_d (cplx_t x, double dr, double di)
/* x = dr + I di */
{
  cplx_Re (x) = dr;
  cplx_Im (x) = di;
}

int
cplx_set_str (cplx_t x, const char *s)
/* set from string as (re , im) */
{
  return sscanf (s, CPLX_INP_FMT, &cplx_Re (x), &cplx_Im (x)) == 2;
}

void
cplx_get_d (double *dr, double *di, const cplx_t x)
/* *dr = re(x), *di = im(x) */
{
  *dr = cplx_Re (x);
  *di = cplx_Im (x);
}

char *
cplx_get_str (char *s, const cplx_t x)
/* output to string as (re , im) */
{
  if (s == NULL && (s = (char *) mps_malloc (DEF_STR_SIZE)) == NULL)
    return NULL;
  sprintf (s, CPLX_OUT_FMT, cplx_Re (x), cplx_Im (x));
  return s;
}

void
cplx_neg (cplx_t rx, const cplx_t x)
/* rx = -x */
{
  cplx_Re (rx) = -cplx_Re (x);
  cplx_Im (rx) = -cplx_Im (x);
}

void
cplx_con (cplx_t rx, const cplx_t x)
/* rx = conj(x) */
{
  cplx_Re (rx) = cplx_Re (x);
  cplx_Im (rx) = -cplx_Im (x);
}

double
cplx_smod (const cplx_t x)
/* returns |x|^2 */
{
  return cplx_Re (x) * cplx_Re (x) + cplx_Im (x) * cplx_Im (x);
}

double
cplx_mod (const cplx_t x)
/* returns |x| */
{
  double d;

  if (fabs (cplx_Re (x)) > fabs (cplx_Im (x)))
    {
      d = cplx_Im (x) / cplx_Re (x);
      return fabs (cplx_Re (x)) * sqrt (1.0 + d * d);
    }
  else if (cplx_Im (x) == 0.0)
    return 0.0;
  d = cplx_Re (x) / cplx_Im (x);
  return fabs (cplx_Im (x)) * sqrt (1.0 + d * d);
}

void
cplx_inv (cplx_t rx, const cplx_t x)
/* rx = 1 / x */
{
  double d1, d2;

  if (fabs (cplx_Re (x)) > fabs (cplx_Im (x)))
    {
      d1 = cplx_Im (x) / cplx_Re (x);
      d2 = 1.0 / (cplx_Re (x) * (1.0 + d1 * d1));
      cplx_Re (rx) = d2;
      cplx_Im (rx) = -d2 * d1;
    }
  else
    {
      d1 = cplx_Re (x) / cplx_Im (x);
      d2 = 1.0 / (cplx_Im (x) * (1.0 + d1 * d1));
      cplx_Im (rx) = -d2;
      cplx_Re (rx) = d2 * d1;
    }
}

void
cplx_sqr (cplx_t rx, const cplx_t x)
/* rx = x * x */
{
  double d;                     /* necessary when rx=x1 or rx=x2 */

  d = cplx_Re (x) * cplx_Re (x) - cplx_Im (x) * cplx_Im (x);
  cplx_Im (rx) = 2.0 * cplx_Re (x) * cplx_Im (x);
  cplx_Re (rx) = d;
}

void
cplx_rot (cplx_t rx, const cplx_t x)
/* rx = I x */
{
  double t;

  t = cplx_Re (x);
  cplx_Re (rx) = -cplx_Im (x);
  cplx_Im (rx) = t;
}

void
cplx_flip (cplx_t rx, const cplx_t x)
/* rx = (Im(x), Re(x)) */
{
  double t;

  t = cplx_Re (x);
  cplx_Re (rx) = cplx_Im (x);
  cplx_Im (rx) = t;
}

/**
 * @brief Add the complex number <code>x1</code> and <code>x2</code> and
 * store the result in <code>rx</code>
 *
 * @param rx The place where the result will be stored.
 * @param x1 The left operand of the addition.
 * @param x2 The right operand of the addition.
 */
void
cplx_add (cplx_t rx, const cplx_t x1, const cplx_t x2)
/* rx = x1 + x2 */
{
  cplx_Re (rx) = cplx_Re (x1) + cplx_Re (x2);
  cplx_Im (rx) = cplx_Im (x1) + cplx_Im (x2);
}

/**
 * @brief Substract from <code>x1</code> the value
 * stored in <code>x2</code> and store the result
 * in <code>rx</code>.
 * 
 * @param rx The place where the result will be stored.
 * @param x1 The left operand of the subtraction.
 * @param x2 The right operand of the subtraction.
 */
void
cplx_sub (cplx_t rx, const cplx_t x1, const cplx_t x2)
/* rx = x1 - x2 */
{
  cplx_Re (rx) = cplx_Re (x1) - cplx_Re (x2);
  cplx_Im (rx) = cplx_Im (x1) - cplx_Im (x2);
}

/**
 * @brief Multiply the value stored in <code>x1</code>
 * and in <code>x2</code> and store the result in
 * <code>rx</code>.
 *
 * @param rx The place where the result wil be stored.
 * @param x1 The left operand of the multiplication.
 * @param x2 The right operando of the multiplication.
 */
void
cplx_mul (cplx_t rx, const cplx_t x1, const cplx_t x2)
/* rx = x1 * x2 */
{
  double d;                     /* necessary when rx=x1 or rx=x2 */

  d = cplx_Re (x1) * cplx_Re (x2) - cplx_Im (x1) * cplx_Im (x2);
  cplx_Im (rx) = cplx_Im (x1) * cplx_Re (x2) + cplx_Re (x1) * cplx_Im (x2);
  cplx_Re (rx) = d;
}

void
cplx_mul_d (cplx_t rx, const cplx_t x, double d)
/* rx = x * d */
{
  cplx_Re (rx) = cplx_Re (x) * d;
  cplx_Im (rx) = cplx_Im (x) * d;
}

/**
 * @brief Divide the complex value in <code>x1</code> by the
 * value in <code>x2</code> and store the results in <code>rx</code>.
 *
 * @param rx The place where the result will be stored.
 * @param x1 The left operand of the division.
 * @param x2 The right operando of the division.
 */
void
cplx_div (cplx_t rx, const cplx_t x1, const cplx_t x2)
  /* rx = x1 / x2 */
{
  cplx_t ctmp;

  cplx_inv (ctmp, x2);
  cplx_mul (rx, x1, ctmp);
}

void
cplx_div_d (cplx_t rx, const cplx_t x, double d)
/* rx = x / d */
{
  cplx_Re (rx) = cplx_Re (x) / d;
  cplx_Im (rx) = cplx_Im (x) / d;
}

void
cplx_pow_si (cplx_t rx, const cplx_t x, register signed long int i)
/* rx = x ^ i , i integer */
{
  cplx_t t;

  cplx_Move (t, x);
  cplx_Move (rx, cplx_one);

  if (i < 0)
    {
      cplx_inv_eq (t);
      i = -i;
    }
  while (i)
    {
      if (i & 1)
        cplx_mul_eq (rx, t);
      cplx_sqr_eq (t);
      i >>= 1;                  /* divide i by 2 */
    }
}

void
cplx_swap (cplx_t x1, cplx_t x2)
/* x1 <-> x2 */
{
  cplx_t t;

  cplx_Move (t, x1);
  cplx_Move (x1, x2);
  cplx_Move (x2, t);
}

void
cplx_neg_eq (cplx_t x)
/* x = -x */
{
  cplx_Re (x) = -cplx_Re (x);
  cplx_Im (x) = -cplx_Im (x);
}

void
cplx_con_eq (cplx_t x)
/* x = conj(x) */
{
  cplx_Im (x) = -cplx_Im (x);
}

void
cplx_inv_eq (cplx_t x)
  /* rx = 1 / x */
{
  double d1, d2;

  if (fabs (cplx_Re (x)) > fabs (cplx_Im (x)))
    {
      d1 = cplx_Im (x) / cplx_Re (x);
      if (DBL_MAX / (1.0 + d1 * d1) < fabs (cplx_Re (x)))       /*#G 27/4/98 */
        d2 = 0.0;               /*#G 27/4/98 */
      else                      /*#G 27/4/98 */
        d2 = 1.0 / (cplx_Re (x) * (1.0 + d1 * d1));
      cplx_Re (x) = d2;
      cplx_Im (x) = -d2 * d1;
    }
  else
    {
      d1 = cplx_Re (x) / cplx_Im (x);
      if (DBL_MAX / (1.0 + d1 * d1) < fabs (cplx_Re (x)))       /*#G 27/4/98 */
        d2 = 0.0;               /*#G 27/4/98 */
      else                      /*#G 27/4/98 */
        d2 = 1.0 / (cplx_Im (x) * (1.0 + d1 * d1));
      cplx_Im (x) = -d2;
      cplx_Re (x) = d2 * d1;
    }
}

void
cplx_sqr_eq (cplx_t x)
/* x = x * x */
{
  double d;

  d = cplx_Re (x) * cplx_Re (x) - cplx_Im (x) * cplx_Im (x);
  cplx_Im (x) *= 2.0 * cplx_Re (x);
  cplx_Re (x) = d;
}

void
cplx_rot_eq (cplx_t x)
/* x = I x */
{
  double d;

  d = cplx_Re (x);
  cplx_Re (x) = -cplx_Im (x);
  cplx_Im (x) = d;
}

void
cplx_flip_eq (cplx_t x)
/* x = (Im(x), Re(x)) */
{
  double d;

  d = cplx_Re (x);
  cplx_Re (x) = cplx_Im (x);
  cplx_Im (x) = d;
}

void
cplx_add_eq (cplx_t rx, const cplx_t x)
/* rx = rx + x */
{
  cplx_Re (rx) += cplx_Re (x);
  cplx_Im (rx) += cplx_Im (x);
}

void
cplx_sub_eq (cplx_t rx, const cplx_t x)
/* rx = rx - x */
{
  cplx_Re (rx) -= cplx_Re (x);
  cplx_Im (rx) -= cplx_Im (x);
}

void
cplx_mul_eq (cplx_t rx, const cplx_t x)
/* rx = rx * x */
{
  double d;

  d = cplx_Re (rx) * cplx_Re (x) - cplx_Im (rx) * cplx_Im (x);
  cplx_Im (rx) = cplx_Im (rx) * cplx_Re (x) + cplx_Re (rx) * cplx_Im (x);
  cplx_Re (rx) = d;
}

void
cplx_mul_eq_d (cplx_t x, double d)
/* x = x * d */
{
  cplx_Re (x) *= d;
  cplx_Im (x) *= d;
}

void
cplx_div_eq (cplx_t rx, const cplx_t x)
  /* rx = x1 / x2 */
{
  cplx_t ctmp;

  cplx_inv (ctmp, x);
  cplx_mul_eq (rx, ctmp);
}

void
cplx_div_eq_d (cplx_t x, double d)
/* x = x / d */
{
  cplx_Re (x) /= d;
  cplx_Im (x) /= d;
}

void
cplx_pow_eq_si (cplx_t x, register signed long int i)
/* x = x ^ i , i integer */
{
  cplx_t t;

  cplx_Move (t, x);
  cplx_Move (x, cplx_one);

  if (i < 0)
    {
      cplx_inv_eq (t);
      i = -i;
    }
  while (i)
    {
      if (i & 1)
        cplx_mul_eq (x, t);
      cplx_sqr_eq (t);
      i >>= 1;                  /* divide i by 2 */
    }
}

/*------------  relational ops.  -------------------------*/

int
cplx_eq_zero (const cplx_t x)
/* x == 0 */
{
  return cplx_Re (x) == 0.0 && cplx_Im (x) == 0.0;
}

int
cplx_eq (const cplx_t x1, const cplx_t x2)
/* x1 == x2 */
{
  return cplx_Re (x1) == cplx_Re (x2) && cplx_Im (x1) == cplx_Im (x2);
}

int
cplx_ne (const cplx_t x1, const cplx_t x2)
/* x1 != x2 */
{
  return cplx_Re (x1) != cplx_Re (x2) || cplx_Im (x1) != cplx_Im (x2);
}

/*------------  vector functions  ------------------------*/

void
cplx_vinit (cplx_t v[], long size)
{
  long i;

  for (i = 0; i < size; i++)
    cplx_Move (v[i], cplx_zero);
}

#else

const cplx_t cplx_zero = { 0.0 };
const cplx_t cplx_one = { 1.0 };
const cplx_t cplx_i = { 1.0I };

#endif

int
cplx_check_fpe (cplx_t c)
/* Check if the components are NaN or Inf */
{
  int fp = 0;
  if (isnan (cplx_Re (c)))
    fp += 1;
  if (isnan (cplx_Im (c)))
    fp += (1 << 1);
  if (isinf (cplx_Re (c)))
    fp += (1 << 2);
  if (isinf (cplx_Im(c)))
    fp += (1 << 3);
  return fp;
}

/*------------  I/O functions  ---------------------------*/

int
cplx_out_str_u (FILE * f, const cplx_t x)
/* output to file as re im */
{
  if (f == NULL)
    f = stdout;
  return fprintf (f, CPLX_OUT_UFMT, cplx_Re (x), cplx_Im (x));
}

int
cplx_out_str (FILE * f, const cplx_t x)
/* output to file as (re, im) */
{
  if (f == NULL)
    f = stdout;
  return fprintf (f, CPLX_OUT_FMT, cplx_Re (x), cplx_Im (x));
}

int
cplx_inp_str_u (cplx_t x, FILE * f)
/* input from file as (re, im) */
{
  double real, imag;
  if (f == NULL)
    f = stdin;
  int ret = fscanf (f, CPLX_INP_UFMT, &real, &imag);
  cplx_set_d (x, real, imag);
  return ret;
}

int
cplx_inp_str (cplx_t x, FILE * f)
/* input from file as (re, im) */
{
  double real, imag;
  if (f == NULL)
    f = stdin;
  int ret = fscanf (f, CPLX_INP_FMT, &real, &imag);
  cplx_set_d (x, real, imag);
  return ret;
}

/***********************************************************
**              functions for rdpe_t                      **
***********************************************************/

/* normalize: mant(e) in [1/2, 1) */
#define rdpe_Norm(E) \
{ \
  int i; \
  rdpe_Mnt(E) = frexp(rdpe_Mnt(E), &i); \
  if (rdpe_Mnt(E) == 0.0) rdpe_Esp(E) = 0L; \
  else rdpe_Esp(E) += i; \
}

/* constants */
const rdpe_t rdpe_zero = { {0.0, 0L} };
const rdpe_t rdpe_one = { {0.5, 1L} };
const rdpe_t RDPE_MAX = { {0.5, LONG_MAX} };
const rdpe_t RDPE_MIN = { {0.5, LONG_MIN} };
const rdpe_t rdpe_maxd = { {0.5, DBL_MAX_EXP} };
const rdpe_t rdpe_mind = { {0.5, DBL_MIN_EXP} };
const rdpe_t RDPE_BIG = { {0.5, LONG_MAX >> 10 } };

void
rdpe_d (rdpe_t temp_rdpe, double d)
/* return d as a rdpe_t */
{
  rdpe_Mnt (temp_rdpe) = d;
  rdpe_Esp (temp_rdpe) = 0L;
  rdpe_Norm (temp_rdpe);
}

void
rdpe_2dl (rdpe_t temp_rdpe, double d, long l)
/* return d*2^l as a rdpe_t */
{
  rdpe_Mnt (temp_rdpe) = d;
  rdpe_Esp (temp_rdpe) = l;
  rdpe_Norm (temp_rdpe);;
}

void
rdpe_clear (rdpe_t e)
/* e = 0 */
{
  rdpe_Move (e, rdpe_zero);
}

void
rdpe_set (rdpe_t re, const rdpe_t e)
/* re = e */
{
  rdpe_Move (re, e);
}

void
rdpe_set_dl (rdpe_t e, double d, long int l)
/* e = d*10^l */
{
  double x;

  if (d == 0.0)
    {
      rdpe_Move (e, rdpe_zero);
      return;
    }
  else if (d > 0.0)
    {
      rdpe_Mnt (e) = log (d) / LOG_2 + l * LOG2_10;
      rdpe_Mnt (e) = pow (2.0, modf (rdpe_Mnt (e), &x));
    }
  else
    {                           /* d < 0 */
      rdpe_Mnt (e) = log (-d) / LOG_2 + l * LOG2_10;
      rdpe_Mnt (e) = -pow (2.0, modf (rdpe_Mnt (e), &x));
    }
  rdpe_Esp (e) = (long int) x;
  rdpe_Norm (e);
}

void
rdpe_set_2dl (rdpe_t e, double d, long int l)
/* e = d*2^l */
{
  rdpe_Mnt (e) = d;
  rdpe_Esp (e) = l;
  rdpe_Norm (e);
}

void
rdpe_set_d (rdpe_t e, double d)
/* e = d  */
{
  rdpe_Mnt (e) = d;
  rdpe_Esp (e) = 0L;
  rdpe_Norm (e);
}

int
rdpe_set_str (rdpe_t e, const char *s)
/* input from string as mant x exp) */
{
  if (sscanf (s, RDPE_INP_FMT, &rdpe_Mnt (e), &rdpe_Esp (e)) != 2)
    return 0;
  rdpe_set_dl (e, rdpe_Mnt (e), rdpe_Esp (e));
  return 1;
}

void
rdpe_get_dl (double *d, long int *l, const rdpe_t e)
/* returns mantissa and exponent of e */
{
  double x;

  if (rdpe_Mnt (e) == 0.0)
    {
      *d = 0.0;
      *l = 0L;
    }
  else if (rdpe_Mnt (e) > 0)
    {
      *d = log10 (rdpe_Mnt (e)) + rdpe_Esp (e) * LOG10_2;
      *d = pow (10.0, modf (*d, &x));
      *l = (long int) x;
    }
  else
    {                           /* rdpe_Mnt(e) > 0 */
      *d = log10 (-rdpe_Mnt (e)) + rdpe_Esp (e) * LOG10_2;
      *d = -pow (10.0, modf (*d, &x));
      *l = (long int) x;
    }
}

void
rdpe_get_2dl (double *d, long int *l, const rdpe_t e)
/* returns mantissa and exponent of e in base 2 */
{
  *d = rdpe_Mnt (e);
  *l = rdpe_Esp (e);
}

double
rdpe_get_d (const rdpe_t e)
{
  return ldexp (rdpe_Mnt (e), (int) rdpe_Esp (e));
}

char *
rdpe_get_str (char *s, const rdpe_t e)
/* output to string */
{
  double d;
  long int l;

  if (s == NULL && (s = (char *) mps_malloc (DEF_STR_SIZE)) == NULL)
    return NULL;
  rdpe_get_dl (&d, &l, e);
  sprintf (s, RDPE_OUT_FMT, d, l);
  return s;
}

void
rdpe_neg (rdpe_t re, const rdpe_t e)
/* re = -e */
{
  rdpe_Mnt (re) = -rdpe_Mnt (e);
  rdpe_Esp (re) = rdpe_Esp (e);
}

void
rdpe_abs (rdpe_t re, const rdpe_t e)
/* re = |e| */
{
  rdpe_Mnt (re) = (rdpe_Mnt (e) > 0) ? rdpe_Mnt (e) : -rdpe_Mnt (e);
  rdpe_Esp (re) = rdpe_Esp (e);
}

void
rdpe_inv (rdpe_t re, const rdpe_t e)
/* re = 1 / e */
{
  rdpe_Mnt (re) = 1.0 / rdpe_Mnt (e);
  rdpe_Esp (re) = -rdpe_Esp (e);
  rdpe_Norm (re);
}

void
rdpe_sqr (rdpe_t re, const rdpe_t e)
/* re = e * e */
{
  rdpe_Mnt (re) = rdpe_Mnt (e) * rdpe_Mnt (e);
  rdpe_Esp (re) = rdpe_Esp (e) + rdpe_Esp (e);
  rdpe_Norm (re);
}

void
rdpe_sqrt (rdpe_t re, const rdpe_t e)
/* re = e^(1/2) */
{
  if (rdpe_Esp (e) & 1)
    {
      rdpe_Mnt (re) = sqrt (rdpe_Mnt (e) / 2.0);
      rdpe_Esp (re) = (rdpe_Esp (e) + 1) / 2;
    }
  else
    {
      rdpe_Mnt (re) = sqrt (rdpe_Mnt (e));
      rdpe_Esp (re) = rdpe_Esp (e) / 2;
    }
  rdpe_Norm (re);
}

double
rdpe_log (const rdpe_t e)
/* returns log(e) */
{
  return log (rdpe_Mnt (e)) + rdpe_Esp (e) * LOG_2;
}

double
rdpe_log10 (const rdpe_t e)
/* returns log10(e) */
{
  return (log (rdpe_Mnt (e)) + rdpe_Esp (e) * LOG_2) / LOG_10;
}

void
rdpe_exp (rdpe_t re, const rdpe_t e)
/* re = E^(e) */
{
  long int i;

  i = rdpe_Esp (e);
  rdpe_set_2dl (re, exp (rdpe_Mnt (e)), 0L);

  if (i >= 0)
    while (i > 0)
      {
        rdpe_sqr_eq (re);
        i--;
      }
  else
    while (i < 0)
      {
        rdpe_sqrt_eq (re);
        i++;
      }
}

void
rdpe_mul (rdpe_t re, const rdpe_t e1, const rdpe_t e2)
/* re = e1 * e2 */
{
  if (rdpe_Esp (e1) >= 0 && (rdpe_Esp (e2) >= LONG_MAX - rdpe_Esp (e1)))
    {
      rdpe_set (re, RDPE_MAX);
      return;
    }
  if (rdpe_Esp (e1) <= 0 && (rdpe_Esp (e2) <= LONG_MIN - rdpe_Esp (e1)))
    {
      rdpe_set (re, RDPE_MAX);
      return;
    }
  rdpe_Mnt (re) = rdpe_Mnt (e1) * rdpe_Mnt (e2);
  rdpe_Esp (re) = rdpe_Esp (e1) + rdpe_Esp (e2);
  rdpe_Norm (re);
}

void
rdpe_mul_d (rdpe_t re, const rdpe_t e, double d)
/* re = e * d */
{
  int esp;
  frexp (d, &esp);
  if (rdpe_Esp (e) >= 0 && (esp >= LONG_MAX - rdpe_Esp (e)))
    {
      rdpe_set (re, RDPE_MAX);
      return;
    }
  if (rdpe_Esp (e) <= 0 && (esp <= LONG_MIN - rdpe_Esp (e)))
    {
      rdpe_set (re, RDPE_MAX);
      return;
    }
  rdpe_Mnt (re) = rdpe_Mnt (e) * d;
  rdpe_Esp (re) = rdpe_Esp (e);
  rdpe_Norm (re);
}

void
rdpe_mul_2exp (rdpe_t re, const rdpe_t e, unsigned long int i)
/* re = e * 2^i */
{
  rdpe_Mnt (re) = rdpe_Mnt (e);
  rdpe_Esp (re) = rdpe_Esp (e) + i;
}

void
rdpe_div (rdpe_t re, const rdpe_t e1, const rdpe_t e2)
/* re = e1 / e2 */
{
  rdpe_Mnt (re) = rdpe_Mnt (e1) / rdpe_Mnt (e2);
  rdpe_Esp (re) = rdpe_Esp (e1) - rdpe_Esp (e2);
  rdpe_Norm (re);
}

void
rdpe_div_d (rdpe_t re, const rdpe_t e, double d)
/* re = e / d */
{
  rdpe_Mnt (re) = rdpe_Mnt (e) / d;
  rdpe_Esp (re) = rdpe_Esp (e);
  rdpe_Norm (re);
}

void
rdpe_div_2exp (rdpe_t re, const rdpe_t e, unsigned long int i)
/* re = e / 2^i */
{
  rdpe_Mnt (re) = rdpe_Mnt (e);
  rdpe_Esp (re) = rdpe_Esp (e) - i;
}

void
rdpe_add (rdpe_t re, const rdpe_t e1, const rdpe_t e2)
/* re = e1 + e2 */
{
  long delta;

  /* Check for overflows */
  if (rdpe_Mnt (e1) > 0 && rdpe_Mnt (e2) > 0 && rdpe_Esp (e1) == LONG_MAX && rdpe_Esp (e2) == LONG_MAX)
    {
      rdpe_set (re, RDPE_MAX);
      return;
    }
  if (rdpe_Mnt (e1) < 0 && rdpe_Mnt (e2) < 0 && rdpe_Esp (e1) == LONG_MAX && rdpe_Esp (e2) == LONG_MAX)
    {
      rdpe_set_dl (re, -0.5, LONG_MAX);
      return;
    }

  if (rdpe_Mnt (e2) == 0.0)
    {
      rdpe_Move (re, e1);
      return;
    }
  if (rdpe_Mnt (e1) == 0.0)
    {
      rdpe_Move (re, e2);
      return;
    }
  delta = rdpe_Esp (e1) - rdpe_Esp (e2);

  if (delta > NBT)
    rdpe_Move (re, e1);
  else if (delta < -NBT)
    rdpe_Move (re, e2);
  else if (delta == 0)
    {
      rdpe_Mnt (re) = rdpe_Mnt (e1) + rdpe_Mnt (e2);
      rdpe_Esp (re) = rdpe_Esp (e1);
      rdpe_Norm (re);
    }
  else if (delta > 0)
    {
      rdpe_Mnt (re) = rdpe_Mnt (e1) + ldexp (rdpe_Mnt (e2), (int) -delta);
      rdpe_Esp (re) = rdpe_Esp (e1);
      rdpe_Norm (re);
    }
  else
    {                           /* delta < 0 */
      rdpe_Mnt (re) = ldexp (rdpe_Mnt (e1), (int) delta) + rdpe_Mnt (e2);
      rdpe_Esp (re) = rdpe_Esp (e2);
      rdpe_Norm (re);
    }
}

void
rdpe_sub (rdpe_t re, const rdpe_t e1, const rdpe_t e2)
/* re = e1 - e2 */
{
  long delta;

  if (rdpe_Mnt (e2) == 0.0)
    {
      rdpe_Move (re, e1);
      return;
    }
  if (rdpe_Mnt (e1) == 0.0)
    {
      rdpe_Mnt (re) = -rdpe_Mnt (e2);
      rdpe_Esp (re) = rdpe_Esp (e2);
      return;
    }

  delta = rdpe_Esp (e1) - rdpe_Esp (e2);

  if (delta > NBT)
    rdpe_Move (re, e1);
  else if (delta < -NBT)
    {
      rdpe_Mnt (re) = -rdpe_Mnt (e2);
      rdpe_Esp (re) = rdpe_Esp (e2);
    }
  else if (delta == 0)
    {
      rdpe_Mnt (re) = rdpe_Mnt (e1) - rdpe_Mnt (e2);
      rdpe_Esp (re) = rdpe_Esp (e1);
      rdpe_Norm (re);
    }
  else if (delta > 0)
    {
      rdpe_Mnt (re) = rdpe_Mnt (e1) - ldexp (rdpe_Mnt (e2), (int) -delta);
      rdpe_Esp (re) = rdpe_Esp (e1);
      rdpe_Norm (re);
    }
  else
    {                           /* delta < 0 */
      rdpe_Mnt (re) = ldexp (rdpe_Mnt (e1), (int) delta) - rdpe_Mnt (e2);
      rdpe_Esp (re) = rdpe_Esp (e2);
      rdpe_Norm (re);
    }
}

void
rdpe_add_d (rdpe_t re, const rdpe_t e, double d)
/* re = e + d */
{
  rdpe_t t;

  rdpe_set_d (t, d);
  rdpe_add (re, e, t);
}

void
rdpe_sub_d (rdpe_t re, const rdpe_t e, double d)
/* re = e - d */
{
  rdpe_t t;

  rdpe_set_d (t, d);
  rdpe_sub (re, e, t);
}

void
rdpe_pow_d (rdpe_t re, const rdpe_t e, double d)
/* re = e ^ d */
{
  double a, i, f;

  a = d * (log (rdpe_Mnt (e)) / LOG_2 + rdpe_Esp (e));
  f = modf (a, &i);
  rdpe_set_2dl (re, exp (f * LOG_2), (long) i);
}

void
rdpe_pow_d2 (rdpe_t re, const rdpe_t e, double d)
/* re = e ^ d, d double */
{
  double a, i;

  a = pow (rdpe_Mnt (e), d);
  a *= pow (2.0, modf (rdpe_Esp (e) * d, &i));
  rdpe_set_2dl (re, a, (long) i);
}

void
rdpe_pow_si (rdpe_t re, const rdpe_t e, register signed long int i)
/* re = e ^ i, i integer */
{
  rdpe_t t;

  rdpe_Move (t, e);
  rdpe_Move (re, rdpe_one);

  if (i < 0)
    {
      rdpe_inv (t, t);
      i = -i;
    }
  while (i)
    {
      if (i & 1)
        rdpe_mul_eq (re, t);
      rdpe_sqr_eq (t);
      i >>= 1;                  /* divide i by 2 */
    }
}

void
rdpe_swap (rdpe_t e1, rdpe_t e2)
/* e1 <-> e2 */
{
  rdpe_t t;

  rdpe_Move (t, e1);
  rdpe_Move (e1, e2);
  rdpe_Move (e2, t);
}

void
rdpe_neg_eq (rdpe_t e)
/* e = -e */
{
  rdpe_Mnt (e) = -rdpe_Mnt (e);
  rdpe_Esp (e) = rdpe_Esp (e);
}

void
rdpe_abs_eq (rdpe_t e)
/* e = |e| */
{
  rdpe_Mnt (e) = (rdpe_Mnt (e) > 0) ? rdpe_Mnt (e) : -rdpe_Mnt (e);
}

void
rdpe_inv_eq (rdpe_t e)
/* e = 1 / e */
{
  rdpe_Mnt (e) = 1.0 / rdpe_Mnt (e);
  rdpe_Esp (e) = -rdpe_Esp (e);
  rdpe_Norm (e);
}

void
rdpe_sqr_eq (rdpe_t e)
/* e = e * e */
{
  rdpe_Mnt (e) *= rdpe_Mnt (e);
  rdpe_Esp (e) <<= 1;           /* multiply by 2 */
  rdpe_Norm (e);
}

void
rdpe_sqrt_eq (rdpe_t e)
/* e = e^(1/2) */
{
  if (rdpe_Esp (e) & 1)
    {                           /* odd test */
      rdpe_Mnt (e) = sqrt (rdpe_Mnt (e) / 2.0);
      rdpe_Esp (e) = (rdpe_Esp (e) + 1) / 2;
    }
  else
    {
      rdpe_Mnt (e) = sqrt (rdpe_Mnt (e));
      rdpe_Esp (e) /= 2;
    }
  rdpe_Norm (e);
}

void
rdpe_exp_eq (rdpe_t e)
/* re = E^(e) */
{
  long int i;

  i = rdpe_Esp (e);
  rdpe_set_2dl (e, exp (rdpe_Mnt (e)), 0L);

  if (i >= 0)
    while (i > 0)
      {
        rdpe_sqr_eq (e);
        i--;
      }
  else
    while (i < 0)
      {
        rdpe_sqrt_eq (e);
        i++;
      }
}

void
rdpe_mul_eq (rdpe_t re, const rdpe_t e)
/* re = re * e */
{;
  if (rdpe_Esp (re) >= 0 && (rdpe_Esp (e) >= LONG_MAX - rdpe_Esp (re)))
    {
      rdpe_set (re, RDPE_MAX);
      return;
    }
  if (rdpe_Esp (re) <= 0 && (rdpe_Esp (e) <= LONG_MIN - rdpe_Esp (re)))
    {
      rdpe_set (re, RDPE_MAX);
      return;
    }
  rdpe_Mnt (re) *= rdpe_Mnt (e);
  rdpe_Esp (re) += rdpe_Esp (e);
  rdpe_Norm (re);
}

void
rdpe_mul_eq_d (rdpe_t e, double d)
/* e = e * d */
{
  int esp;
  frexp (d, &esp);
  if (rdpe_Esp (e) >= 0 && (esp >= LONG_MAX - rdpe_Esp (e)))
    {
      rdpe_set (e, RDPE_MAX);
      return;
    }
  if (rdpe_Esp (e) <= 0 && (esp <= LONG_MIN - rdpe_Esp (e)))
    {
      rdpe_set (e, RDPE_MAX);
      return;
    }
  rdpe_Mnt (e) *= d;
  rdpe_Norm (e);
}

void
rdpe_mul_eq_2exp (rdpe_t e, unsigned long int i)
/* e = e * 2^i */
{
  rdpe_Esp (e) += i;
}

void
rdpe_div_eq (rdpe_t re, const rdpe_t e)
/* re = re / e */
{
  rdpe_Mnt (re) /= rdpe_Mnt (e);
  rdpe_Esp (re) -= rdpe_Esp (e);
  rdpe_Norm (re);
}

void
rdpe_div_eq_d (rdpe_t e, double d)
/* e = e / d */
{
  rdpe_Mnt (e) /= d;
  rdpe_Norm (e);
}

void
rdpe_div_eq_2exp (rdpe_t e, unsigned long int i)
/* e = e / 2^i */
{
  rdpe_Esp (e) -= i;
}

void
rdpe_add_eq (rdpe_t re, const rdpe_t e)
/* re = re + e */
{
  long int delta;

  if (rdpe_Mnt (e) == 0.0)
    return;
  if (rdpe_Mnt (re) == 0.0)
    {
      rdpe_Move (re, e);
      return;
    }
  delta = rdpe_Esp (re) - rdpe_Esp (e);

  if (delta > NBT)
    return;
  else if (delta < -NBT)
    rdpe_Move (re, e);
  else if (delta == 0)
    {
      rdpe_Mnt (re) += rdpe_Mnt (e);
      rdpe_Norm (re);
    }
  else if (delta > 0)
    {
      rdpe_Mnt (re) += ldexp (rdpe_Mnt (e), (int) -delta);
      rdpe_Norm (re);
    }
  else
    {                           /* delta < 0 */
      rdpe_Mnt (re) = ldexp (rdpe_Mnt (re), (int) delta) + rdpe_Mnt (e);
      rdpe_Esp (re) = rdpe_Esp (e);
      rdpe_Norm (re);
    }
}

void
rdpe_sub_eq (rdpe_t re, const rdpe_t e)
/* re = re - e */
{
  long int delta;

  if (rdpe_Mnt (e) == 0.0)
    return;
  if (rdpe_Mnt (re) == 0.0)
    {
      rdpe_Mnt (re) = -rdpe_Mnt (e);
      rdpe_Esp (re) = rdpe_Esp (e);
      return;
    }
  delta = rdpe_Esp (re) - rdpe_Esp (e);

  if (delta > NBT)
    return;
  else if (delta < -NBT)
    {
      rdpe_Mnt (re) = -rdpe_Mnt (e);
      rdpe_Esp (re) = rdpe_Esp (e);
    }
  else if (delta == 0)
    {
      rdpe_Mnt (re) -= rdpe_Mnt (e);
      rdpe_Norm (re);
    }
  else if (delta > 0)
    {
      rdpe_Mnt (re) -= ldexp (rdpe_Mnt (e), (int) -delta);
      rdpe_Norm (re);
    }
  else
    {                           /* delta < 0 */
      rdpe_Mnt (re) = ldexp (rdpe_Mnt (re), (int) delta) - rdpe_Mnt (e);
      rdpe_Esp (re) = rdpe_Esp (e);
      rdpe_Norm (re);
    }
}

void
rdpe_add_eq_d (rdpe_t e, double d)
/* re = e + d */
{
  rdpe_t t;

  rdpe_set_d (t, d);
  rdpe_add_eq (e, t);
}

void
rdpe_sub_eq_d (rdpe_t e, double d)
/* re = e - d */
{
  rdpe_t t;

  rdpe_set_d (t, d);
  rdpe_sub_eq (e, t);
}

void
rdpe_pow_eq_d (rdpe_t e, double d)
/* re = e ^ d */
{
  double a, i, f;

  a = d * (log (rdpe_Mnt (e)) / LOG_2 + rdpe_Esp (e));
  f = modf (a, &i);
  rdpe_set_2dl (e, exp (f * LOG_2), (long) i);
}

void
rdpe_pow_eq_si (rdpe_t e, register signed long int i)
/* e = e ^ i, i integer */
{
  rdpe_t t;

  rdpe_Move (t, e);
  rdpe_Move (e, rdpe_one);

  if (i < 0)
    {
      rdpe_inv (t, t);
      i = -i;
    }
  while (i)
    {
      if (i & 1)
        rdpe_mul_eq (e, t);
      rdpe_sqr_eq (t);
      i >>= 1;                  /* divide i by 2 */
    }
}

void
rdpe_fac_ui (rdpe_t e, register unsigned long int n)
/* e = n! */
{
  rdpe_Move (e, rdpe_one);
  while (n > 1)
    {
      rdpe_mul_eq_d (e, (double) n);
      n--;
    }
}

/*------------  relational ops.  -------------------------*/

int
rdpe_cmp (const rdpe_t e1, const rdpe_t e2)
/* e1 <!> e2 */
{
  rdpe_t t;

  rdpe_sub (t, e1, e2);
  if (rdpe_Mnt (t) > 0.0)
    return 1;
  if (rdpe_Mnt (t) < 0.0)
    return -1;
  return 0;
}

int
rdpe_sgn (const rdpe_t e)
/* sign(e) */
{
  if (rdpe_Mnt (e) > 0.0)
    return 1;
  if (rdpe_Mnt (e) < 0.0)
    return -1;
  return 0;
}

int
rdpe_eq_zero (const rdpe_t e)
/* e == 0 */
{
  return rdpe_Mnt (e) == 0.0 && rdpe_Esp (e) == 0;
}

int
rdpe_eq (const rdpe_t e1, const rdpe_t e2)
/* e1 == e2 */
{
  return rdpe_Mnt (e1) == rdpe_Mnt (e2) && rdpe_Esp (e1) == rdpe_Esp (e2);
}

int
rdpe_ne (const rdpe_t e1, const rdpe_t e2)
/* e1 != e2 */
{
  return rdpe_Mnt (e1) != rdpe_Mnt (e2) || rdpe_Esp (e1) != rdpe_Esp (e2);
}

int
rdpe_lt (const rdpe_t e1, const rdpe_t e2)
/* e1 < e2 */
{
  rdpe_t t;

  if (rdpe_Mnt (e1) > 0 && rdpe_Mnt (e2) < 0)
      return 1;
  if (rdpe_Mnt (e1) < 0 && rdpe_Mnt (e1) > 0)
      return 0;

  /* This check works only if the numbers are non zero */
  if (rdpe_Mnt (e1) != 0 && rdpe_Mnt (e2) != 0)
    {
      if (rdpe_Esp (e1) > rdpe_Esp (e2))
          return 0;
      if (rdpe_Esp (e2) > rdpe_Esp (e1)) 
          return 1; 
    }

  rdpe_sub (t, e1, e2);
  return rdpe_Mnt (t) < 0.0;
}

int
rdpe_le (const rdpe_t e1, const rdpe_t e2)
/* e1 <= e2 */
{
  rdpe_t t;

  if (rdpe_Mnt (e1) > 0 && rdpe_Mnt (e2) < 0)
    return 1;
  if (rdpe_Mnt (e1) < 0 && rdpe_Mnt (e1) > 0)
    return 0;

  /* This check works only if the numbers are non zero */
  if (rdpe_Mnt (e1) != 0 && rdpe_Mnt (e2) != 0)
    {
      if (rdpe_Esp (e1) > rdpe_Esp (e2)) 
        return 0;
      if (rdpe_Esp (e2) > rdpe_Esp (e1)) 
        return 1;
    }

  rdpe_sub (t, e1, e2);
  return rdpe_Mnt (t) <= 0.0;
}

int
rdpe_gt (const rdpe_t e1, const rdpe_t e2)
/* e1 > e2 */
{
  rdpe_t t;

  if (rdpe_Mnt (e1) > 0 && rdpe_Mnt (e2) < 0)
    return 0;
  if (rdpe_Mnt (e1) < 0 && rdpe_Mnt (e1) > 0)
    return 1;

  /* This check works only if both DPE are non zero */
  if (rdpe_Mnt (e1) != 0 && rdpe_Mnt (e2) != 0)
    {       
      if (rdpe_Esp (e1) > rdpe_Esp (e2)) 
        return 1; 
      if (rdpe_Esp (e2) > rdpe_Esp (e1)) 
        return 0; 
    }

  rdpe_sub (t, e1, e2);
  return rdpe_Mnt (t) > 0.0;
}

int
rdpe_ge (const rdpe_t e1, const rdpe_t e2)
/* e1 >= e2 */
{
  rdpe_t t;

  if (rdpe_Mnt (e1) > 0 && rdpe_Mnt (e2) < 0)
    return 0;
  if (rdpe_Mnt (e1) < 0 && rdpe_Mnt (e1) > 0)
    return 1;

  /* This check works only if both DPE are non zero */
  if (rdpe_Mnt (e1) != 0 && rdpe_Mnt (e2) != 0)
    {       
      if (rdpe_Esp (e1) > rdpe_Esp (e2)) 
        return 1; 
      if (rdpe_Esp (e2) > rdpe_Esp (e1)) 
        return 0; 
    }

  rdpe_sub (t, e1, e2);
  return rdpe_Mnt (t) >= 0.0;
}

/*------------  I/O functions  ---------------------------*/

int
rdpe_out_str_u (FILE * f, const rdpe_t e)
/* output as mantissa e exponent, base 10 */
{
  double d;
  long int l;

  if (f == NULL)
    f = stdout;
  rdpe_get_dl (&d, &l, e);
  return fprintf (f, RDPE_OUT_UFMT, d, l);
}

int
rdpe_out_str (FILE * f, const rdpe_t e)
/* output as mantissa x exponent, base 10 */
{
  double d;
  long int l;

  if (f == NULL)
    f = stdout;
  rdpe_get_dl (&d, &l, e);
  return fprintf (f, RDPE_OUT_FMT, d, l);
}

int
rdpe_inp_str_u (rdpe_t e, FILE * f)
/* input from file as mant e exp */
{
  double d;
  long int l;

  if (f == NULL)
    f = stdin;
  if (fscanf (f, RDPE_INP_UFMT, &d, &l) != 2)
    return 0;
  rdpe_set_dl (e, d, l);
  return 1;
}

int
rdpe_inp_str (rdpe_t e, FILE * f)
/* input from file as mant x exp */
{
  double d;
  long int l;

  if (f == NULL)
    f = stdin;
  if (fscanf (f, RDPE_INP_FMT, &d, &l) != 2)
    return 0;
  rdpe_set_dl (e, d, l);
  return 1;
}

int
rdpe_inp_str_flex (rdpe_t e, FILE * f)
/* More flexible input for rdpe */
{
  double d;
  long int l = 0;

  if (f == NULL)
    f = stdin;

  if (fscanf (f, RDPE_INP_FMT, &d, &l) < 1)
    return 0;
  rdpe_set_dl (e, d, l);
  return 1;
}

int
rdpe_inp_sstr_flex (rdpe_t e, char *f)
/* More flexible input for rdpe */
{
  double d;
  long int l = 0;

  if (sscanf (f, RDPE_INP_FMT, &d, &l) < 1)
    return 0;
  rdpe_set_dl (e, d, l);
  return 1;
}

/*------------  vector functions  ------------------------*/

void
rdpe_vinit (rdpe_t v[], long size)
{
  long i;

  for (i = 0; i < size; i++)
    rdpe_Move (v[i], rdpe_zero);
}

void 
gdpe_add (gdpe_t res, gdpe_t g1, gdpe_t g2)
{
  rdpe_add (gdpe_Val (res), gdpe_Val (g1), gdpe_Val (g2));
  rdpe_add (gdpe_Eps (res), gdpe_Eps (g1), gdpe_Eps (g2));
  gdpe_update_rel_from_abs (res);
}

void 
gdpe_sub (gdpe_t res, gdpe_t g1, gdpe_t g2)
{
  rdpe_sub (gdpe_Val (res), gdpe_Val (g1), gdpe_Val (g2));
  rdpe_add (gdpe_Eps (res), gdpe_Eps (g1), gdpe_Eps (g2));
  gdpe_update_rel_from_abs (res);
}

void 
gdpe_mul (gdpe_t res, gdpe_t g1, gdpe_t g2)
{
  rdpe_mul (gdpe_Val (res), gdpe_Val (g1), gdpe_Val (g2));
  if (gdpe_eq_zero (g1) || gdpe_eq_zero (g2))
    {
      rdpe_set (gdpe_Eps (res), rdpe_zero);
      return;
    }

  rdpe_add (gdpe_Rel (res), gdpe_Rel (g1), gdpe_Rel (g2));
  gdpe_update_abs_from_rel (res);
}

void 
gdpe_div (gdpe_t res, gdpe_t g1, gdpe_t g2)
{
  rdpe_div (gdpe_Val (res), gdpe_Val (g1), gdpe_Val (g2));
  if (gdpe_eq_zero (g1))
    {
      rdpe_set (gdpe_Eps (res), rdpe_zero);
      return;
    }

  rdpe_add (gdpe_Rel (res), gdpe_Rel (g1), gdpe_Rel (g2));
  gdpe_update_abs_from_rel (res);
}


/***********************************************************
**              functions for cdpe_t                      **
***********************************************************/

/* normalize both parts: mantissa in [1/2, 1) */

#define cdpe_Norm(C)  rdpe_Norm(cdpe_Re(C)); rdpe_Norm(cdpe_Im(C));

/* base constants */
const cdpe_t cdpe_zero = { {{{0.0, 0L}}, {{0.0, 0L}}} };
const cdpe_t cdpe_one = { {{{0.5, 1L}}, {{0.0, 0L}}} };
const cdpe_t cdpe_i = { {{{0.0, 0L}}, {{0.5, 1L}}} };

void
cdpe_d (cdpe_t temp_cdpe, double r, double i)
/* return (r, i) as a cdpe_t */
{
  rdpe_Mnt (cdpe_Re (temp_cdpe)) = r;
  rdpe_Esp (cdpe_Re (temp_cdpe)) = 0L;
  rdpe_Mnt (cdpe_Im (temp_cdpe)) = i;
  rdpe_Esp (cdpe_Im (temp_cdpe)) = 0L;
  cdpe_Norm (temp_cdpe);
}

void
cdpe_x (cdpe_t temp_cdpe, const cplx_t x)
/* return x as a cdpe_t */
{
  rdpe_Mnt (cdpe_Re (temp_cdpe)) = cplx_Re (x);
  rdpe_Esp (cdpe_Re (temp_cdpe)) = 0L;
  rdpe_Mnt (cdpe_Im (temp_cdpe)) = cplx_Im (x);
  rdpe_Esp (cdpe_Im (temp_cdpe)) = 0L;
  cdpe_Norm (temp_cdpe);
}

void
cdpe_e (cdpe_t temp_cdpe, const rdpe_t er, const rdpe_t ei)
/* return (er, ei) */
{
  rdpe_Move (cdpe_Re (temp_cdpe), er);
  rdpe_Move (cdpe_Im (temp_cdpe), ei);
}

void
cdpe_2dl (cdpe_t temp_cdpe, double dr, long lr, double di, long li)
/* return (dr*2^lr, di*2^li) as a cdpe_t */
{
  rdpe_Mnt (cdpe_Re (temp_cdpe)) = dr;
  rdpe_Esp (cdpe_Re (temp_cdpe)) = lr;
  rdpe_Mnt (cdpe_Im (temp_cdpe)) = di;
  rdpe_Esp (cdpe_Im (temp_cdpe)) = li;
  cdpe_Norm (temp_cdpe);
}

void
cdpe_init (cdpe_t c)
/* c = 0 + I 0 */
{
  cdpe_Move (c, cdpe_zero);
}

void
cdpe_set (cdpe_t rc, const cdpe_t c)
/* rc = c */
{
  cdpe_Move (rc, c);
}

void
cdpe_set_x (cdpe_t c, const cplx_t x)
/* c = (cdpe_t) x */
{
  cdpe_Move (c, cdpe_zero);
  rdpe_Mnt (cdpe_Re (c)) = cplx_Re (x);
  rdpe_Mnt (cdpe_Im (c)) = cplx_Im (x);
  cdpe_Norm (c);
}

void
cdpe_set_e (cdpe_t c, const rdpe_t er, const rdpe_t ei)
/* c = er + I ei */
{
  rdpe_Move (cdpe_Re (c), er);
  rdpe_Move (cdpe_Im (c), ei);
}

void
cdpe_set_d (cdpe_t c, double dr, double di)
/* c = dr + I di */
{
  cdpe_set (c, cdpe_zero);
  rdpe_Mnt (cdpe_Re (c)) = dr;
  rdpe_Mnt (cdpe_Im (c)) = di;
  cdpe_Norm (c);
}

void
cdpe_set_dl (cdpe_t c, double dr, long int lr, double di, long int li)
/* c = dr*10^lr + I di*10^li */
{
  rdpe_set_dl (cdpe_Re (c), dr, lr);
  rdpe_set_dl (cdpe_Im (c), di, li);
}

void
cdpe_set_2dl (cdpe_t c, double dr, long int lr, double di, long int li)
/* c = dr*2^lr + I di*2^li */
{
  rdpe_Mnt (cdpe_Re (c)) = dr;
  rdpe_Esp (cdpe_Re (c)) = lr;
  rdpe_Mnt (cdpe_Im (c)) = di;
  rdpe_Esp (cdpe_Im (c)) = li;
  cdpe_Norm (c);
}

int
cdpe_set_str (cdpe_t c, const char *s)
/* set from string as (re , im) */
{
  if (sscanf
      (s, CDPE_INP_FMT, &rdpe_Mnt (cdpe_Re (c)), &rdpe_Esp (cdpe_Re (c)),
       &rdpe_Mnt (cdpe_Im (c)), &rdpe_Esp (cdpe_Im (c))) != 4)
    return 0;
  rdpe_set_dl (cdpe_Re (c), rdpe_Mnt (cdpe_Re (c)), rdpe_Esp (cdpe_Re (c)));
  rdpe_set_dl (cdpe_Im (c), rdpe_Mnt (cdpe_Im (c)), rdpe_Esp (cdpe_Im (c)));
  return 1;
}

void
cdpe_get_e (rdpe_t er, rdpe_t ei, const cdpe_t c)
/* er = re(c), ei = im(c) */
{
  rdpe_Move (er, cdpe_Re (c));
  rdpe_Move (ei, cdpe_Im (c));
}

void
cdpe_get_x (cplx_t x, const cdpe_t c)
/* e = im(c) */
{
  cplx_set_d (x, 
              ldexp (rdpe_Mnt (cdpe_Re (c)), (int) rdpe_Esp (cdpe_Re (c))),
              ldexp (rdpe_Mnt (cdpe_Im (c)), (int) rdpe_Esp (cdpe_Im (c))));
}

void
cdpe_get_d (double *dr, double *di, const cdpe_t c)
/* *dr = re(c), *di = im(c) */
{
  *dr = ldexp (rdpe_Mnt (cdpe_Re (c)), (int) rdpe_Esp (cdpe_Re (c)));
  *di = ldexp (rdpe_Mnt (cdpe_Im (c)), (int) rdpe_Esp (cdpe_Im (c)));
}

char *
cdpe_get_str (char *s, const cdpe_t c)
/* output to string as (re , im) */
{
  double dr, di;
  long int lr, li;

  if (s == NULL && (s = (char *) mps_malloc (DEF_STR_SIZE)) == NULL)
    return NULL;
  rdpe_get_dl (&dr, &lr, cdpe_Re (c));
  rdpe_get_dl (&di, &li, cdpe_Im (c));
  sprintf (s, CDPE_OUT_FMT, dr, lr, di, li);
  return s;
}

void
cdpe_neg (cdpe_t rc, const cdpe_t c)
/* rc = -c */
{
  cdpe_Move (rc, c);
  rdpe_Mnt (cdpe_Re (rc)) = -rdpe_Mnt (cdpe_Re (rc));
  rdpe_Mnt (cdpe_Im (rc)) = -rdpe_Mnt (cdpe_Im (rc));
}

void
cdpe_con (cdpe_t rc, const cdpe_t c)
/* rc = conj(c) */
{
  cdpe_Move (rc, c);
  rdpe_Mnt (cdpe_Im (rc)) = -rdpe_Mnt (cdpe_Im (rc));
}

void
cdpe_smod (rdpe_t e, const cdpe_t c)
/* e = |c|^2 */
{
  rdpe_t t;

  rdpe_sqr (e, cdpe_Re (c));
  rdpe_sqr (t, cdpe_Im (c));
  rdpe_add_eq (e, t);
}

void
cdpe_mod (rdpe_t e, const cdpe_t c)
/* e = |c| */
{
  rdpe_t t;

  rdpe_sqr (e, cdpe_Re (c));
  rdpe_sqr (t, cdpe_Im (c));
  rdpe_add_eq (e, t);
  rdpe_sqrt_eq (e);
}

void
cdpe_inv (cdpe_t rc, const cdpe_t c)
/* rc = 1 / c */
{
  rdpe_t e;

  cdpe_smod (e, c);
  rdpe_inv_eq (e);
  cdpe_Move (rc, c);
  rdpe_Mnt (cdpe_Im (rc)) = -rdpe_Mnt (cdpe_Im (rc));
  rdpe_mul_eq (cdpe_Re (rc), e);
  rdpe_mul_eq (cdpe_Im (rc), e);
}

void
cdpe_sqr (cdpe_t rc, const cdpe_t c)
/* rc = c * c */
{
  rdpe_t e1, e2;

  rdpe_mul (e1, cdpe_Re (c), cdpe_Re (c));
  rdpe_mul (e2, cdpe_Im (c), cdpe_Im (c));
  rdpe_mul (cdpe_Im (rc), cdpe_Im (c), cdpe_Re (c));
  rdpe_Esp (cdpe_Im (rc)) += 1; /* multiply by 2 */
  rdpe_sub (cdpe_Re (rc), e1, e2);
}

void
cdpe_rot (cdpe_t rc, const cdpe_t c)
/* rc = I c */
{
  rdpe_t e;

  rdpe_Move (e, cdpe_Re (c));
  rdpe_Move (cdpe_Re (rc), cdpe_Im (c));
  rdpe_Move (cdpe_Im (rc), e);
  rdpe_Mnt (cdpe_Re (rc)) = -rdpe_Mnt (cdpe_Re (rc));
}

void
cdpe_flip (cdpe_t rc, const cdpe_t c)
/* rc = (Im(c), Re(c)) */
{
  rdpe_t e;

  rdpe_Move (e, cdpe_Re (c));
  rdpe_Move (cdpe_Re (rc), cdpe_Im (c));
  rdpe_Move (cdpe_Im (rc), e);
}

void
cdpe_add (cdpe_t rc, const cdpe_t c1, const cdpe_t c2)
/* rc = c1 + c2 */
{
  rdpe_add (cdpe_Re (rc), cdpe_Re (c1), cdpe_Re (c2));
  rdpe_add (cdpe_Im (rc), cdpe_Im (c1), cdpe_Im (c2));
}

void
cdpe_sub (cdpe_t rc, const cdpe_t c1, const cdpe_t c2)
/* rc = c1 - c2 */
{
  rdpe_sub (cdpe_Re (rc), cdpe_Re (c1), cdpe_Re (c2));
  rdpe_sub (cdpe_Im (rc), cdpe_Im (c1), cdpe_Im (c2));
}

void
cdpe_mul (cdpe_t rc, const cdpe_t c1, const cdpe_t c2)
/* rc = c1 * c2 */
{
  rdpe_t e1, e2, e3;

  rdpe_mul (e1, cdpe_Re (c1), cdpe_Re (c2));
  rdpe_mul (e2, cdpe_Im (c1), cdpe_Im (c2));
  rdpe_sub (e3, e1, e2);        /* needed when rc=c1 or rc=c2 */
  rdpe_mul (e1, cdpe_Im (c1), cdpe_Re (c2));
  rdpe_mul (e2, cdpe_Re (c1), cdpe_Im (c2));
  rdpe_Move (cdpe_Re (rc), e3);
  rdpe_add (cdpe_Im (rc), e1, e2);
}

void
cdpe_mul_e (cdpe_t rc, const cdpe_t c, const rdpe_t e)
/* rc = c * e */
{
  rdpe_Mnt (cdpe_Re (rc)) = rdpe_Mnt (cdpe_Re (c)) * rdpe_Mnt (e);
  rdpe_Esp (cdpe_Re (rc)) = rdpe_Esp (cdpe_Re (c)) + rdpe_Esp (e);
  rdpe_Mnt (cdpe_Im (rc)) = rdpe_Mnt (cdpe_Im (c)) * rdpe_Mnt (e);
  rdpe_Esp (cdpe_Im (rc)) = rdpe_Esp (cdpe_Im (c)) + rdpe_Esp (e);
  cdpe_Norm (rc);
}

void
cdpe_mul_x (cdpe_t rc, const cdpe_t c, const cplx_t x)
/* rc = c * x */
{
  rdpe_t e1, e2, e3;

  rdpe_mul_d (e1, cdpe_Re (c), cplx_Re (x));
  rdpe_mul_d (e2, cdpe_Im (c), cplx_Im (x));
  rdpe_sub (e3, e1, e2);        /* needed when rc=c1 or rc=c2 */
  rdpe_mul_d (e1, cdpe_Im (c), cplx_Re (x));
  rdpe_mul_d (e2, cdpe_Re (c), cplx_Im (x));
  rdpe_Move (cdpe_Re (rc), e3);
  rdpe_add (cdpe_Im (rc), e1, e2);
}

void
cdpe_mul_d (cdpe_t rc, const cdpe_t c, double d)
/* rc = c * d */
{
  cdpe_Move (rc, c);
  rdpe_Mnt (cdpe_Re (rc)) *= d;
  rdpe_Mnt (cdpe_Im (rc)) *= d;
  cdpe_Norm (rc);
}

void
cdpe_mul_2exp (cdpe_t rc, const cdpe_t c, unsigned long int i)
/* rc = c * 2^i */
{
  cdpe_Move (rc, c);
  rdpe_Esp (cdpe_Re (rc)) += i;
  rdpe_Esp (cdpe_Im (rc)) += i;
}

void
cdpe_div (cdpe_t rc, const cdpe_t c1, const cdpe_t c2)
/* rc = c1 / c2 */
{
  cdpe_t t;                     /* needed when rc=c1 or rc=c2 */
  rdpe_t e1, e2, e3;

  cdpe_smod (e1, c2);
  cdpe_div_e (t, c2, e1);
  rdpe_Mnt (cdpe_Im (t)) = -rdpe_Mnt (cdpe_Im (t));
  rdpe_mul (e1, cdpe_Re (c1), cdpe_Re (t));
  rdpe_mul (e2, cdpe_Im (c1), cdpe_Im (t));
  rdpe_sub (e3, e1, e2);
  rdpe_mul (e1, cdpe_Im (c1), cdpe_Re (t));
  rdpe_mul (e2, cdpe_Re (c1), cdpe_Im (t));
  rdpe_Move (cdpe_Re (rc), e3);
  rdpe_add (cdpe_Im (rc), e1, e2);
}

void
cdpe_div_e (cdpe_t rc, const cdpe_t c, const rdpe_t e)
/* rc = c / e */
{
  rdpe_Mnt (cdpe_Re (rc)) = rdpe_Mnt (cdpe_Re (c)) / rdpe_Mnt (e);
  rdpe_Esp (cdpe_Re (rc)) = rdpe_Esp (cdpe_Re (c)) - rdpe_Esp (e);
  rdpe_Mnt (cdpe_Im (rc)) = rdpe_Mnt (cdpe_Im (c)) / rdpe_Mnt (e);
  rdpe_Esp (cdpe_Im (rc)) = rdpe_Esp (cdpe_Im (c)) - rdpe_Esp (e);
  cdpe_Norm (rc);
}

void
cdpe_div_d (cdpe_t rc, const cdpe_t c, double d)
/* rc = c / d */
{
  cdpe_Move (rc, c);
  rdpe_Mnt (cdpe_Re (rc)) /= d;
  rdpe_Mnt (cdpe_Im (rc)) /= d;
  cdpe_Norm (rc);
}

void
cdpe_div_2exp (cdpe_t rc, const cdpe_t c, unsigned long int i)
/* rc = c / 2^i */
{
  cdpe_Move (rc, c);
  rdpe_Esp (cdpe_Re (rc)) -= i;
  rdpe_Esp (cdpe_Im (rc)) -= i;
}

void
cdpe_pow_si (cdpe_t rc, const cdpe_t c, register signed long int i)
/* rc = c^i, i integer */
{
  cdpe_t t;

  cdpe_Move (t, c);
  cdpe_Move (rc, cdpe_one);

  if (i < 0)
    {
      cdpe_inv (t, t);
      i = -i;
    }
  while (i)
    {
      if (i & 1)
        cdpe_mul_eq (rc, t);
      cdpe_sqr_eq (t);
      i >>= 1;                  /* divide i by 2 */
    }
}

void
cdpe_swap (cdpe_t c1, cdpe_t c2)
/* c1 <-> c2 */
{
  cdpe_t t;

  cdpe_Move (t, c1);
  cdpe_Move (c1, c2);
  cdpe_Move (c2, t);
}

void
cdpe_neg_eq (cdpe_t c)
/* c = -c */
{
  rdpe_Mnt (cdpe_Re (c)) = -rdpe_Mnt (cdpe_Re (c));
  rdpe_Mnt (cdpe_Im (c)) = -rdpe_Mnt (cdpe_Im (c));
}

void
cdpe_con_eq (cdpe_t c)
/* c = conj(c) */
{
  rdpe_Mnt (cdpe_Im (c)) = -rdpe_Mnt (cdpe_Im (c));
}

void
cdpe_inv_eq (cdpe_t c)
/* c = 1 / c */
{
  rdpe_t e;

  cdpe_smod (e, c);
  rdpe_inv_eq (e);
  rdpe_Mnt (cdpe_Im (c)) = -rdpe_Mnt (cdpe_Im (c));
  rdpe_mul_eq (cdpe_Re (c), e);
  rdpe_mul_eq (cdpe_Im (c), e);
}

void
cdpe_sqr_eq (cdpe_t c)
/* c = c * c */
{
  rdpe_t e1, e2;

  rdpe_sqr (e1, cdpe_Re (c));
  rdpe_sqr (e2, cdpe_Im (c));
  rdpe_mul_eq (cdpe_Im (c), cdpe_Re (c));
  rdpe_Esp (cdpe_Im (c)) += 1;  /* multiply by 2 */
  rdpe_sub (cdpe_Re (c), e1, e2);
}

void
cdpe_rot_eq (cdpe_t c)
/* c = I c */
{
  rdpe_t e;

  rdpe_Move (e, cdpe_Re (c));
  rdpe_Move (cdpe_Re (c), cdpe_Im (c));
  rdpe_Move (cdpe_Im (c), e);
  rdpe_Mnt (cdpe_Re (c)) = -rdpe_Mnt (cdpe_Re (c));
}

void
cdpe_flip_eq (cdpe_t c)
/* c = (Im(c), Re(c)) */
{
  rdpe_t e;

  rdpe_Move (e, cdpe_Re (c));
  rdpe_Move (cdpe_Re (c), cdpe_Im (c));
  rdpe_Move (cdpe_Im (c), e);
}

void
cdpe_add_eq (cdpe_t rc, const cdpe_t c)
/* rc = rc + c2 */
{
  rdpe_add_eq (cdpe_Re (rc), cdpe_Re (c));
  rdpe_add_eq (cdpe_Im (rc), cdpe_Im (c));
}

void
cdpe_sub_eq (cdpe_t rc, const cdpe_t c)
/* rc = rc - c2 */
{
  rdpe_sub_eq (cdpe_Re (rc), cdpe_Re (c));
  rdpe_sub_eq (cdpe_Im (rc), cdpe_Im (c));
}

void
cdpe_mul_eq (cdpe_t rc, const cdpe_t c)
/* rc = rc * c */
{
  rdpe_t e1, e2, e3;

  rdpe_mul (e1, cdpe_Re (rc), cdpe_Re (c));
  rdpe_mul (e2, cdpe_Im (rc), cdpe_Im (c));
  rdpe_sub (e3, e1, e2);
  rdpe_mul (e1, cdpe_Im (rc), cdpe_Re (c));
  rdpe_mul (e2, cdpe_Re (rc), cdpe_Im (c));
  rdpe_Move (cdpe_Re (rc), e3);
  rdpe_add (cdpe_Im (rc), e1, e2);
}

void
cdpe_mul_eq_e (cdpe_t c, const rdpe_t e)
/* c = c * e */
{
  rdpe_Mnt (cdpe_Re (c)) *= rdpe_Mnt (e);
  rdpe_Esp (cdpe_Re (c)) += rdpe_Esp (e);
  rdpe_Mnt (cdpe_Im (c)) *= rdpe_Mnt (e);
  rdpe_Esp (cdpe_Im (c)) += rdpe_Esp (e);
  cdpe_Norm (c);
}

void
cdpe_mul_eq_x (cdpe_t c, const cplx_t x)
/* c *= x */
{
  rdpe_t e1, e2, e3;

  rdpe_mul_d (e1, cdpe_Re (c), cplx_Re (x));
  rdpe_mul_d (e2, cdpe_Im (c), cplx_Im (x));
  rdpe_sub (e3, e1, e2);        /* needed when rc=c1 or rc=c2 */
  rdpe_mul_d (e1, cdpe_Im (c), cplx_Re (x));
  rdpe_mul_d (e2, cdpe_Re (c), cplx_Im (x));
  rdpe_Move (cdpe_Re (c), e3);
  rdpe_add (cdpe_Im (c), e1, e2);
}

void
cdpe_mul_eq_d (cdpe_t c, double d)
/* c = c * d */
{
  rdpe_Mnt (cdpe_Re (c)) *= d;
  rdpe_Mnt (cdpe_Im (c)) *= d;
  cdpe_Norm (c);
}

void
cdpe_mul_eq_2exp (cdpe_t c, unsigned long int i)
/* c = c * 2^i */
{
  rdpe_Esp (cdpe_Re (c)) += i;
  rdpe_Esp (cdpe_Im (c)) += i;
}

void
cdpe_div_eq (cdpe_t rc, const cdpe_t c)
/* rc = rc / c */
{
  cdpe_t t;
  rdpe_t e1, e2, e3;

  cdpe_smod (e1, c);
  cdpe_div_e (t, c, e1);
  rdpe_Mnt (cdpe_Im (t)) = -rdpe_Mnt (cdpe_Im (t));
  rdpe_mul (e1, cdpe_Re (c), cdpe_Re (t));
  rdpe_mul (e2, cdpe_Im (c), cdpe_Im (t));
  rdpe_sub (e3, e1, e2);
  rdpe_mul (e1, cdpe_Im (c), cdpe_Re (t));
  rdpe_mul (e2, cdpe_Re (c), cdpe_Im (t));
  rdpe_Move (cdpe_Re (rc), e3);
  rdpe_add (cdpe_Im (rc), e1, e2);
}

void
cdpe_div_eq_e (cdpe_t c, const rdpe_t e)
/* c = c / e */
{
  rdpe_Mnt (cdpe_Re (c)) /= rdpe_Mnt (e);
  rdpe_Esp (cdpe_Re (c)) -= rdpe_Esp (e);
  rdpe_Mnt (cdpe_Im (c)) /= rdpe_Mnt (e);
  rdpe_Esp (cdpe_Im (c)) -= rdpe_Esp (e);
  cdpe_Norm (c);
}

void
cdpe_div_eq_d (cdpe_t c, double d)
/* c = c / d */
{
  rdpe_Mnt (cdpe_Re (c)) /= d;
  rdpe_Mnt (cdpe_Im (c)) /= d;
  cdpe_Norm (c);
}

void
cdpe_div_eq_2exp (cdpe_t c, unsigned long int i)
/* c = c / 2^i */
{
  rdpe_Esp (cdpe_Re (c)) -= i;
  rdpe_Esp (cdpe_Im (c)) -= i;
}

void
cdpe_pow_eq_si (cdpe_t c, register signed long int i)
/* rc = c^i, i integer */
{
  cdpe_t t;

  cdpe_Move (t, c);
  cdpe_Move (c, cdpe_one);

  if (i < 0)
    {
      cdpe_inv (t, t);
      i = -i;
    }
  while (i)
    {
      if (i & 1)
        cdpe_mul_eq (c, t);
      cdpe_sqr_eq (t);
      i >>= 1;                  /* divide i by 2 */
    }
}

/*------------  relational ops.  -------------------------*/

int
cdpe_eq_zero (const cdpe_t c)
/* c == 0 */
{
  return rdpe_Mnt (cdpe_Re (c)) == 0.0 && rdpe_Esp (cdpe_Re (c)) == 0 &&
    rdpe_Mnt (cdpe_Im (c)) == 0.0 && rdpe_Esp (cdpe_Im (c)) == 0;
}

int
cdpe_eq (const cdpe_t c1, const cdpe_t c2)
/* c1 == c2 */
{
  return rdpe_Mnt (cdpe_Re (c1)) == rdpe_Mnt (cdpe_Re (c2)) &&
    rdpe_Esp (cdpe_Re (c1)) == rdpe_Esp (cdpe_Re (c2)) &&
    rdpe_Mnt (cdpe_Im (c1)) == rdpe_Mnt (cdpe_Im (c2)) &&
    rdpe_Esp (cdpe_Im (c1)) == rdpe_Esp (cdpe_Im (c2));
}

int
cdpe_ne (const cdpe_t c1, const cdpe_t c2)
/* c1 != c2 */
{
  return rdpe_Mnt (cdpe_Re (c1)) != rdpe_Mnt (cdpe_Re (c2)) ||
    rdpe_Esp (cdpe_Re (c1)) != rdpe_Esp (cdpe_Re (c2)) ||
    rdpe_Mnt (cdpe_Im (c1)) != rdpe_Mnt (cdpe_Im (c2)) ||
    rdpe_Esp (cdpe_Im (c1)) != rdpe_Esp (cdpe_Im (c2));
}

/*------------  I/O functions  ---------------------------*/

int
cdpe_out_str_u (FILE * f, const cdpe_t c)
/* output as mant e exp  mant e exp  */
{
  if (f == NULL)
    f = stdout;
  if (rdpe_out_str_u (f, cdpe_Re (c)) < 0)
    return -1;
  if (fprintf (f, " ") < 0)
    return -1;
  if (rdpe_out_str_u (f, cdpe_Im (c)) < 0)
    return -1;
  return 0;
}

int
cdpe_out_str (FILE * f, const cdpe_t c)
/* output as (mant x exp , mant x exp)  */
{
  if (f == NULL)
    f = stdout;
  if (fputc ('(', f) == EOF)
    return -1;
  if (rdpe_out_str (f, cdpe_Re (c)) < 0)
    return -1;
  if (fprintf (f, ", ") < 0)
    return -1;
  if (rdpe_out_str (f, cdpe_Im (c)) < 0)
    return -1;
  return fputc (')', f);
}

int
cdpe_inp_str_u (cdpe_t c, FILE * f)
/* input from file as mant e exp  mant e exp */
{
  double dr, di;
  long int lr, li;

  if (f == NULL)
    f = stdin;
  if (fscanf (f, CDPE_INP_UFMT, &dr, &lr, &di, &li) != 4)
    return 0;
  rdpe_set_dl (cdpe_Re (c), dr, lr);
  rdpe_set_dl (cdpe_Im (c), di, li);
  return 1;
}

int
cdpe_inp_str (cdpe_t c, FILE * f)
/* input from file as (mant x exp , mant x exp) */
{
  double dr, di;
  long int lr, li;

  if (f == NULL)
    f = stdin;
  if (fscanf (f, CDPE_INP_FMT, &dr, &lr, &di, &li) != 4)
    return 0;
  rdpe_set_dl (cdpe_Re (c), dr, lr);
  rdpe_set_dl (cdpe_Im (c), di, li);
  return 1;
}

/*------------  vector functions  ------------------------*/

void
cdpe_vinit (cdpe_t v[], long size)
{
  long i;

  for (i = 0; i < size; i++)
    cdpe_Move (v[i], cdpe_zero);
}

/***********************************************************
**                                                        **
***********************************************************/
