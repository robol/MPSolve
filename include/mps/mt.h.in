/***********************************************************
**      More Types library for extended C arithmetic      **
**                      Version 1.1                       **
**                                                        **
**             Written by Giuseppe Fiorentino             **
**                 (fiorent@dm.unipi.it)                  **
**                                                        **
** (C) 1997, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

/**
 * @file
 * @brief Library with extended types in C.
 * 
 * It contains complex types based on double and
 * DPE types (Double Plus Exponent) in real and
 * complex versions.
 */

#ifndef __MT_H__
#define __MT_H__

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

/***********************************************************
**              cplx_t type                               **
***********************************************************/

#ifdef MPS_USE_BUILTIN_COMPLEX

/* ----   definition  ---- */
  typedef struct
  {
    double r, i;                /* re, im */
  } __cplx_struct;

  typedef __cplx_struct cplx_t[1];
  typedef const __cplx_struct *cplx_cp;

/* macros for fields access */

#define cplx_Val(X)       (*(X))
#define cplx_Re(X)        ((X)->r)
#define cplx_Im(X)        ((X)->i)
#define cplx_Addr(X)      ((__cplx_struct *) (X))
#define cplx_Move(X1, X2) (*(X1) = *(X2))

/* base constants */
  extern const cplx_t cplx_zero;        /* complex zero  (0, 0)    */
  extern const cplx_t cplx_one; /* complex one   (1, 0)    */
  extern const cplx_t cplx_i;   /* imag. unit    (0, 1)    */

/* built constants */
  void cplx_d (cplx_t temp_cplx, double r, double i);

/* initializers */
#define cplx_init(X) cplx_clear(X)
  void cplx_clear (cplx_t x);
  void cplx_set (cplx_t rx, const cplx_t x);
  void cplx_set_d (cplx_t x, double dr, double di);
  int cplx_set_str (cplx_t x, const char *s);

  /* check floating point exceptions */
  int cplx_check_fpe (cplx_t x);

/* conversion */
  void cplx_get_d (double *dr, double *di, const cplx_t x);
  char *cplx_get_str (char *s, const cplx_t x);

/* unary functions */
  void cplx_neg (cplx_t rx, const cplx_t x);
  void cplx_con (cplx_t rx, const cplx_t x);
  double cplx_smod (const cplx_t x);
  double cplx_mod (const cplx_t x);
  void cplx_inv (cplx_t rx, const cplx_t x);
  void cplx_sqr (cplx_t rx, const cplx_t x);
  void cplx_rot (cplx_t rx, const cplx_t x);
  void cplx_flip (cplx_t rx, const cplx_t x);

/* binary functions */
  void cplx_add (cplx_t rx, const cplx_t x1, const cplx_t x2);
  void cplx_sub (cplx_t rx, const cplx_t x1, const cplx_t x2);
  void cplx_mul (cplx_t rx, const cplx_t x1, const cplx_t x2);
  void cplx_div (cplx_t rx, const cplx_t x1, const cplx_t x2);
  void cplx_mul_d (cplx_t rx, const cplx_t x, double d);
  void cplx_div_d (cplx_t rx, const cplx_t x, double d);
  void cplx_pow_si (cplx_t rx, const cplx_t x, signed long int i);

/* misc */
  void cplx_swap (cplx_t x1, cplx_t x2);

/* op= style operators */
  void cplx_neg_eq (cplx_t x);
  void cplx_con_eq (cplx_t x);
  void cplx_inv_eq (cplx_t x);
  void cplx_sqr_eq (cplx_t x);
  void cplx_rot_eq (cplx_t x);
  void cplx_flip_eq (cplx_t x);
  void cplx_add_eq (cplx_t rx, const cplx_t x);
  void cplx_sub_eq (cplx_t rx, const cplx_t x);
  void cplx_mul_eq (cplx_t rx, const cplx_t x);
  void cplx_div_eq (cplx_t rx, const cplx_t x);
  void cplx_mul_eq_d (cplx_t x, double d);
  void cplx_div_eq_d (cplx_t x, double d);
  void cplx_pow_eq_si (cplx_t x, signed long int i);

/* relational operators */
  int cplx_eq_zero (const cplx_t x);
  int cplx_eq (const cplx_t x1, const cplx_t x2);
  int cplx_ne (const cplx_t x1, const cplx_t x2);

/* I/O functions */
  int cplx_out_str_u (FILE * f, const cplx_t x);
  int cplx_out_str (FILE * f, const cplx_t x);
  int cplx_inp_str_u (cplx_t x, FILE * f);
  int cplx_inp_str (cplx_t x, FILE * f);
#define cplx_outln_str_u(F, C) cplx_out_str_u(F, C), fputc('\n', F)
#define cplx_outln_str(F, C)   cplx_out_str(F, C), fputc('\n', F)
#define cplx_out_u(C)          cplx_out_str_u(stdout, C)
#define cplx_out(C)            cplx_out_str(stdout, C)
#define cplx_outln_u(C)        cplx_out_str_u(stdout, C), putchar('\n')
#define cplx_outln(C)          cplx_out_str(stdout, C), putchar('\n')
#define cplx_inp_u(C)          cplx_inp_str_u(C, stdin)
#define cplx_inp(C)            cplx_inp_str(C, stdin)

/* vector functions */
#define cplx_valloc(N)       (cplx_t *) malloc((N) * sizeof(cplx_t))
  void cplx_vinit (cplx_t v[], long size);
/* #define cplx_vclear(V)       free(V) */
/* #define cplx_vclear(V, N)    cplx_vinit(V, N) */
#define cplx_vfree(V)        free(V)

#else

  #include <complex.h>
  typedef _Complex double cplx_t[1];

/* macros for fields access */

#define cplx_Re(X)        (creal(*X))
#define cplx_Im(X)        (cimag(*X))
#define cplx_Addr(X)      ((complex double *) X)

/* base constants */
  extern const cplx_t cplx_zero;        /* complex zero  (0, 0)    */
  extern const cplx_t cplx_one; /* complex one   (1, 0)    */
  extern const cplx_t cplx_i;   /* imag. unit    (0, 1)    */

/* initializers */
#define cplx_set(x,y) (*x = *y)
#define cplx_set_d(x,y,z) (*x = y + 1.0I * z)

  /* check floating point exceptions */
  int cplx_check_fpe (cplx_t x);

/* unary functions */
#define cplx_mod(x) (cabs(*x))
#define cplx_inv(x,y) (*x = 1 / *y)
#define cplx_sqr(x,y) (*x = csqrt(*y))

/* binary functions */
#define cplx_add(x,y,z) (*x = *y + *z)
#define cplx_sub(x,y,z) (*x = *y - *z)
#define cplx_mul(x,y,z) (*x = *y * *z)
#define cplx_div(x,y,z) (*x = *y / *z)
#define cplx_mul_d(x,y,z) (*x = *y * z)
#define cplx_div_d(x,y,z) (*x = *y / z)

#define cplx_add_eq(x,y) (*x += *y)
#define cplx_sub_eq(x,y) (*x -= *y)
#define cplx_mul_eq(x,y) (*x *= *y)
#define cplx_div_eq(x,y) (*x /= *y)
#define cplx_mul_eq_d(x,y) (*x *= y)
#define cplx_div_eq_d(x,y) (*x /= y)
#define cplx_inv_eq(x) (*x = 1.0 / *x)

/* relational operators */
#define cplx_eq_zero(x) (*x == 0.0)
#define cplx_eq(x,y) (*x == *y)
#define cplx_ne(x,y) (*x != éy)

/* I/O functions */
  int cplx_out_str_u (FILE * f, const cplx_t x);
  int cplx_out_str (FILE * f, const cplx_t x);
  int cplx_inp_str_u (cplx_t x, FILE * f);
  int cplx_inp_str (cplx_t x, FILE * f);
#define cplx_outln_str_u(F, C) cplx_out_str_u(F, C), fputc('\n', F)
#define cplx_outln_str(F, C)   cplx_out_str(F, C), fputc('\n', F)
#define cplx_out_u(C)          cplx_out_str_u(stdout, C)
#define cplx_out(C)            cplx_out_str(stdout, C)
#define cplx_outln_u(C)        cplx_out_str_u(stdout, C), putchar('\n')
#define cplx_outln(C)          cplx_out_str(stdout, C), putchar('\n')
#define cplx_inp_u(C)          cplx_inp_str_u(C, stdin)
#define cplx_inp(C)            cplx_inp_str(C, stdin)

/* vector functions */
#define cplx_valloc(N)       (cplx_t *) malloc((N) * sizeof(cplx_t))
/* #define cplx_vclear(V)       free(V) */
/* #define cplx_vclear(V, N)    cplx_vinit(V, N) */
#define cplx_vfree(V)        free(V)

#endif


/***********************************************************
**              rdpe_t type                               ** 
***********************************************************/

/* ----   definition  ---- */
  typedef struct
  {
    double m;                   /* mantissa */
    long e;                     /* exponent */
  } __rdpe_struct;

  typedef __rdpe_struct rdpe_t[1];
  typedef const __rdpe_struct *rdpe_cp;

  /* macros for fields access */
#define rdpe_Val(E)       (*E)
#define rdpe_Mnt(E)       (E -> m)
#define rdpe_Esp(E)       (E -> e)
#define rdpe_Addr(E)      ((__rdpe_struct *) E)
#define rdpe_Move(E1, E2) (*E1 = *E2)

  /* base constants */
  extern const rdpe_t rdpe_zero;        /* zero as rdpe num.       */
  extern const rdpe_t rdpe_one; /* one as rdpe num.        */
  extern const rdpe_t rdpe_maxd;        /* max double as rdpe      */
  extern const rdpe_t rdpe_mind;        /* min pos. double as rdpe */
  extern const rdpe_t RDPE_MAX; /* max rdpe number         */
  extern const rdpe_t RDPE_MIN; /* min pos. rdpe  number   */
  extern const rdpe_t RDPE_BIG;

  /* built constants */
  void rdpe_d (rdpe_t temp_rdpe, double d);
  void rdpe_2dl (rdpe_t temp_rdpe, double d, long l);

  /* assignment functions */
#define rdpe_init(E) rdpe_clear(E)
  void rdpe_clear (rdpe_t e);
  void rdpe_set (rdpe_t re, const rdpe_t e);
  void rdpe_set_d (rdpe_t e, double d);
  void rdpe_set_dl (rdpe_t e, double d, long int l);
  void rdpe_set_2dl (rdpe_t e, double d, long int l);
  int rdpe_set_str (rdpe_t e, const char *s);

  /* conversion */
  double rdpe_get_d (const rdpe_t e);
  void rdpe_get_dl (double *d, long int *l, const rdpe_t e);
  void rdpe_get_2dl (double *d, long int *l, const rdpe_t e);
  char *rdpe_get_str (char *s, const rdpe_t e);

/* unary functions */
  void rdpe_neg (rdpe_t re, const rdpe_t e);
  void rdpe_abs (rdpe_t re, const rdpe_t e);
  void rdpe_inv (rdpe_t re, const rdpe_t e);
  void rdpe_sqr (rdpe_t re, const rdpe_t e);
  void rdpe_sqrt (rdpe_t re, const rdpe_t e);
  double rdpe_log (const rdpe_t e);
  double rdpe_log10 (const rdpe_t e);
  void rdpe_exp (rdpe_t re, const rdpe_t e);

/* binary functions */
  void rdpe_add (rdpe_t re, const rdpe_t e1, const rdpe_t e2);
  void rdpe_sub (rdpe_t re, const rdpe_t e1, const rdpe_t e2);
  void rdpe_mul (rdpe_t re, const rdpe_t e1, const rdpe_t e2);
  void rdpe_div (rdpe_t re, const rdpe_t e1, const rdpe_t e2);
  void rdpe_add_d (rdpe_t re, const rdpe_t e, double d);
  void rdpe_sub_d (rdpe_t re, const rdpe_t e, double d);
  void rdpe_mul_d (rdpe_t re, const rdpe_t e, double d);
  void rdpe_mul_2exp (rdpe_t re, const rdpe_t e, unsigned long int i);
  void rdpe_div_d (rdpe_t re, const rdpe_t e, double d);
  void rdpe_div_2exp (rdpe_t re, const rdpe_t e, unsigned long int i);
  void rdpe_pow_d (rdpe_t re, const rdpe_t e, double d);
  void rdpe_pow_si (rdpe_t re, const rdpe_t e, signed long int i);

/* misc */
  void rdpe_fac_ui (rdpe_t e, unsigned long int n);
  void rdpe_swap (rdpe_t e1, rdpe_t e2);

/* op= style operators */
  void rdpe_neg_eq (rdpe_t e);
  void rdpe_abs_eq (rdpe_t e);
  void rdpe_inv_eq (rdpe_t e);
  void rdpe_sqr_eq (rdpe_t e);
  void rdpe_sqrt_eq (rdpe_t e);
  void rdpe_exp_eq (rdpe_t e);
  void rdpe_add_eq (rdpe_t re, const rdpe_t e);
  void rdpe_sub_eq (rdpe_t re, const rdpe_t e);
  void rdpe_mul_eq (rdpe_t re, const rdpe_t e);
  void rdpe_div_eq (rdpe_t re, const rdpe_t e);
  void rdpe_add_eq_d (rdpe_t e, double d);
  void rdpe_sub_eq_d (rdpe_t e, double d);
  void rdpe_mul_eq_d (rdpe_t e, double d);
  void rdpe_mul_eq_2exp (rdpe_t e, unsigned long int i);
  void rdpe_div_eq_d (rdpe_t e, double d);
  void rdpe_div_eq_2exp (rdpe_t e, unsigned long int i);
  void rdpe_pow_eq_d (rdpe_t e, double d);
  void rdpe_pow_eq_si (rdpe_t e, signed long int i);

/* relational ops */
  int rdpe_cmp (const rdpe_t e1, const rdpe_t e2);
  int rdpe_sgn (const rdpe_t e);
  int rdpe_eq_zero (const rdpe_t e);
  int rdpe_eq (const rdpe_t e1, const rdpe_t e2);
  int rdpe_ne (const rdpe_t e1, const rdpe_t e2);
  int rdpe_lt (const rdpe_t e1, const rdpe_t e2);
  int rdpe_le (const rdpe_t e1, const rdpe_t e2);
  int rdpe_gt (const rdpe_t e1, const rdpe_t e2);
  int rdpe_ge (const rdpe_t e1, const rdpe_t e2);

/* I/O functions */
  int rdpe_out_str_u (FILE * f, const rdpe_t e);
  int rdpe_out_str (FILE * f, const rdpe_t e);
  int rdpe_inp_str_u (rdpe_t e, FILE * f);
  int rdpe_inp_str (rdpe_t e, FILE * f);
  int rdpe_inp_str_flex (rdpe_t e, FILE * f);
#define rdpe_outln_str_u(F, E) rdpe_out_str_u(F, E); fputc('\n', F)
#define rdpe_outln_str(F, E)   rdpe_out_str(F, E); fputc('\n', F)
#define rdpe_out_u(E)          rdpe_out_str_u(stdout, E)
#define rdpe_out(E)            rdpe_out_str(stdout, E)
#define rdpe_outln_u(E)        rdpe_out_str_u(stdout, E); putchar('\n')
#define rdpe_outln(E)          rdpe_out_str(stdout, E); putchar('\n')
#define rdpe_inp_u(E)          rdpe_inp_str_u(e, stdin)
#define rdpe_inp(E)            rdpe_inp_str(e, stdin)

/* vector functions */
#define rdpe_valloc(N)       (rdpe_t *) malloc((N) * sizeof(rdpe_t))
  void rdpe_vinit (rdpe_t v[], long size);
/* #define rdpe_vclear(V)       free(V) */
/* #define rdpe_vclear(V, N)    rdpe_vinit(V, N) */
#define rdpe_vfree(V)        free(V)

/***********************************************************
**              gdpe_t functions                          ** 
***********************************************************/
typedef struct 
{
  rdpe_t r;
  rdpe_t eps;
  rdpe_t rel_eps;
} __gdpe_struct;

  typedef __gdpe_struct gdpe_t[1];

#define gdpe_Eps(g) ((g)->eps)
#define gdpe_Rel(g) ((g)->rel_eps)
#define gdpe_Val(g) ((g)->r)
#define gdpe_update_rel_from_abs(g) (rdpe_div (gdpe_Rel (g), gdpe_Eps (g), gdpe_Val (g)))
#define gdpe_update_abs_from_rel(g) (rdpe_mul (gdpe_Eps (g), gdpe_Rel (g), gdpe_Val (g)))

  /* Tests */
#define gdpe_eq_zero(g) (rdpe_eq_zero (gdpe_Val (g)))

#define gdpe_add_eq (g, g2) (gdpe_add ((g), (g), (g2)))
#define gdpe_sub_eq (g, g2) (gdpe_sub ((g), (g), (g2)))
#define gdpe_mul_eq (g, g2) (gdpe_mul ((g), (g), (g2)))
#define gdpe_div_eq (g, g2) (gdpe_div ((g), (g), (g2)))

  /* Binary functions */
  void gdpe_add (gdpe_t res, gdpe_t g1, gdpe_t g2);
  void gdpe_sub (gdpe_t res, gdpe_t g1, gdpe_t g2);
  void gdpe_mul (gdpe_t res, gdpe_t g1, gdpe_t g2);
  void gdpe_div (gdpe_t res, gdpe_t g1, gdpe_t g2);



/***********************************************************
**              cdpe_t functions                          ** 
***********************************************************/

/* ---- definition ---- */
  typedef struct
  {
    rdpe_t r, i;                /* re, im */
  } __cdpe_struct;

  typedef __cdpe_struct cdpe_t[1];
  typedef const __cdpe_struct *cdpe_cp;

/* macros for fields access */

#define cdpe_Val(C)       (*C)
#define cdpe_Re(C)        (C -> r)
#define cdpe_Im(C)        (C -> i)
#define cdpe_Addr(C)      ((__cdpe_struct *) C)
#define cdpe_Move(C1, C2) (*C1 = *C2)

/* base constants */
  extern const cdpe_t cdpe_zero;        /* cdpe zero     (0, 0)    */
  extern const cdpe_t cdpe_one; /* cdpe one      (1, 0)    */
  extern const cdpe_t cdpe_i;   /* cdpe I        (0, 1)    */

/* built constants */
  void cdpe_d (cdpe_t temp_cdpe, double r, double i);
  void cdpe_x (cdpe_t temp_cdpe, const cplx_t x);
  void cdpe_e (cdpe_t temp_cdpe, const rdpe_t er, const rdpe_t ei);
  void cdpe_2dl (cdpe_t temp_cdpe, double dr, long lr, double di, long li);

/* initializers */
#define cdpe_init(C) cdpe_clear(C)
  void cdpe_clear (cdpe_t c);
  void cdpe_set (cdpe_t rc, const cdpe_t c);
  void cdpe_set_e (cdpe_t c, const rdpe_t er, const rdpe_t ei);
  void cdpe_set_x (cdpe_t c, const cplx_t x);
  void cdpe_set_d (cdpe_t c, double dr, double di);
  void cdpe_set_dl (cdpe_t c, double dr, long int lr, double di, long int li);
  void cdpe_set_2dl (cdpe_t c, double dr, long int lr,
                     double di, long int li);
  int cdpe_set_str (cdpe_t c, const char *s);

/* conversion */
  void cdpe_get_e (rdpe_t er, rdpe_t ei, const cdpe_t c);
  void cdpe_get_x (cplx_t x, const cdpe_t c);
  void cdpe_get_d (double *dr, double *di, const cdpe_t c);
  char *cdpe_get_str (char *s, const cdpe_t c);

/* unary functions */
  void cdpe_neg (cdpe_t rc, const cdpe_t c);
  void cdpe_con (cdpe_t rc, const cdpe_t c);
  void cdpe_smod (rdpe_t e, const cdpe_t c);
  void cdpe_mod (rdpe_t e, const cdpe_t c);
  void cdpe_inv (cdpe_t rc, const cdpe_t c);
  void cdpe_sqr (cdpe_t rc, const cdpe_t c);
  void cdpe_rot (cdpe_t rc, const cdpe_t c);
  void cdpe_flip (cdpe_t rc, const cdpe_t c);

/* binary functions */
  void cdpe_add (cdpe_t rc, const cdpe_t c1, const cdpe_t c2);
  void cdpe_sub (cdpe_t rc, const cdpe_t c1, const cdpe_t c2);
  void cdpe_mul (cdpe_t rc, const cdpe_t c1, const cdpe_t c2);
  void cdpe_div (cdpe_t rc, const cdpe_t c1, const cdpe_t c2);
  void cdpe_mul_e (cdpe_t rc, const cdpe_t c, const rdpe_t e);
  void cdpe_mul_x (cdpe_t rc, const cdpe_t c, const cplx_t x);
  void cdpe_mul_d (cdpe_t rc, const cdpe_t c, double d);
  void cdpe_mul_2exp (cdpe_t rc, const cdpe_t c, unsigned long int i);
  void cdpe_div_e (cdpe_t rc, const cdpe_t c, const rdpe_t e);
  void cdpe_div_d (cdpe_t rc, const cdpe_t c, double d);
  void cdpe_div_2exp (cdpe_t rc, const cdpe_t c, unsigned long int i);

/* misc */
  void cdpe_pow_si (cdpe_t rc, const cdpe_t c, signed long int i);

/* op= style operators */
  void cdpe_neg_eq (cdpe_t c);
  void cdpe_con_eq (cdpe_t c);
  void cdpe_inv_eq (cdpe_t c);
  void cdpe_sqr_eq (cdpe_t c);
  void cdpe_rot_eq (cdpe_t c);
  void cdpe_flip_eq (cdpe_t c);
  void cdpe_add_eq (cdpe_t rc, const cdpe_t c);
  void cdpe_sub_eq (cdpe_t rc, const cdpe_t c);
  void cdpe_mul_eq (cdpe_t rc, const cdpe_t c);
  void cdpe_div_eq (cdpe_t rc, const cdpe_t c);
  void cdpe_mul_eq_e (cdpe_t c, const rdpe_t e);
  void cdpe_mul_eq_x (cdpe_t c, const cplx_t x);
  void cdpe_mul_eq_d (cdpe_t c, double d);
  void cdpe_mul_eq_2exp (cdpe_t c, unsigned long int i);
  void cdpe_div_eq_e (cdpe_t c, const rdpe_t e);
  void cdpe_div_eq_d (cdpe_t c, double d);
  void cdpe_div_eq_2exp (cdpe_t c, unsigned long int i);
  void cdpe_pow_eq_si (cdpe_t c, signed long int i);
  void cdpe_swap (cdpe_t c1, cdpe_t c2);

/* relational ops */
  int cdpe_eq_zero (const cdpe_t c);
  int cdpe_eq (const cdpe_t c1, const cdpe_t c2);
  int cdpe_ne (const cdpe_t c1, const cdpe_t c2);

/* I/O functions */
  int cdpe_out_str_u (FILE * f, const cdpe_t c);
  int cdpe_out_str (FILE * f, const cdpe_t c);
  int cdpe_inp_str_u (cdpe_t c, FILE * f);
  int cdpe_inp_str (cdpe_t c, FILE * f);
#define cdpe_outln_str_u(F, C) cdpe_out_str_u(F, C); fputc('\n', F)
#define cdpe_outln_str(F, C)   cdpe_out_str(F, C); fputc('\n', F)
#define cdpe_out_u(C)          cdpe_out_str_u(stdout, C)
#define cdpe_out(C)            cdpe_out_str(stdout, C)
#define cdpe_outln_u(C)        cdpe_out_str_u(stdout, C); putchar('\n')
#define cdpe_outln(C)          cdpe_out_str(stdout, C); putchar('\n')
#define cdpe_inp_u(C)          cdpe_inp_str_u(C, stdin)
#define cdpe_inp(C)            cdpe_inp_str(C, stdin)

/* vector functions */
#define cdpe_valloc(N)       (cdpe_t *) malloc((N) * sizeof(cdpe_t))
  void cdpe_vinit (cdpe_t v[], long size);
/* #define cdpe_vclear(V)       free(C) */
/* #define cdpe_vclear(V, N)    cdpe_vinit(V, N) */
#define cdpe_vfree(V)        free(V)

/*
 * End of extern "C" {
 *   ...
 * }
 */
#ifdef __cplusplus
}
#endif

#endif

/***********************************************************
**                                                        **
***********************************************************/
