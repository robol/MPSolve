/***********************************************************
**          Multi-Precision complex type for GMP          **
**                      Version 1.1                       **
**                                                        **
**             Written by Giuseppe Fiorentino             **
**                 (fiorent@dm.unipi.it)                  **
***********************************************************/

/**
 * @file
 * @brief Multiprecision complex type, based on mpf multiprecision
 * floating point type.
 */

#ifndef __MPC_H__
#define __MPC_H__

#include <mps/mt.h>
#include <stdio.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/***********************************************************
**              definition of types                       **
***********************************************************/

  typedef struct
  {
    mpf_t r, i;
  } __mpc_struct;

  typedef __mpc_struct mpc_t[1];

/***********************************************************
**            macros for fields access                    **
***********************************************************/

/*
 * #define mpc_Val(C)             (*(C))
 */
#define mpc_Re(C)               ((C)->r)
#define mpc_Im(C)               ((C)->i)
#define mpc_Addr(C)             ((__mpc_struct *) C)
#define mpc_Move(C1, C2)        (*(C1) = *(C2))

/***********************************************************
**            mpc_t functions                             **
***********************************************************/

/* constructors/destructors */
  void mpc_init (mpc_t c);
  void mpc_init2 (mpc_t c, unsigned long int prec);
  void mpc_clear (mpc_t c);

  void mpc_set_prec (mpc_t c, unsigned long int prec);
  unsigned long int mpc_get_prec (mpc_t c);
  void mpc_set_prec_raw (mpc_t c, unsigned long int prec);

/* initializers */
  void mpc_set (mpc_t rc, mpc_t c);
  void mpc_set_ui (mpc_t c, unsigned long int ir, unsigned long int ii);
  void mpc_set_si (mpc_t c, signed long int ir, signed long int ii);
  void mpc_set_d (mpc_t c, double dr, double di);
  void mpc_set_z (mpc_t c, mpz_t zr, mpz_t zi);
  void mpc_set_q (mpc_t c, mpq_t qr, mpq_t qi);
  void mpc_set_f (mpc_t c, mpf_t fr, mpf_t fi);
  int mpc_set_str (mpc_t c, char *sr, char *si, int base);

  void mpc_init_set (mpc_t rc, mpc_t c);
  void mpc_init_set_ui (mpc_t c, unsigned long int ir, unsigned long int ii);
  void mpc_init_set_si (mpc_t c, signed long int ir, signed long int ii);
  void mpc_init_set_d (mpc_t c, double dr, double di);
  void mpc_init_set_f (mpc_t c, mpf_t fr, mpf_t fi);
  int mpc_init_set_str (mpc_t c, char *sr, char *si, int base);

/* unary functions */
  void mpc_neg (mpc_t rc, mpc_t c);
  void mpc_smod (mpf_t f, mpc_t c);
  void mpc_rmod (rdpe_t r, mpc_t c);
  void mpc_mod (mpf_t f, mpc_t c);
  void mpc_con (mpc_t rc, mpc_t c);
  void mpc_inv (mpc_t rc, mpc_t c);
  void mpc_inv2 (mpc_t rc, mpc_t c);
  void mpc_sqr (mpc_t rc, mpc_t c);
  void mpc_rot (mpc_t rc, mpc_t c);
  void mpc_flip (mpc_t rc, mpc_t c);

/* binary functions */
  void mpc_add (mpc_t rc, mpc_t c1, mpc_t c2);
  void mpc_add_f (mpc_t rc, mpc_t c, mpf_t f);
  void mpc_add_ui (mpc_t rc, mpc_t c, unsigned long int r,
                   unsigned long int i);
  void mpc_sub (mpc_t rc, mpc_t c1, mpc_t c2);
  void mpc_sub_f (mpc_t rc, mpc_t c, mpf_t f);
  void mpc_f_sub (mpc_t rc, mpf_t f, mpc_t c);
  void mpc_sub_ui (mpc_t rc, mpc_t c, unsigned long int r,
                   unsigned long int i);
  void mpc_ui_sub (mpc_t rc, unsigned long int r, unsigned long int i,
                   mpc_t c);
  void mpc_mul (mpc_t rc, mpc_t c1, mpc_t c2);
  void mpc_mul_f (mpc_t rc, mpc_t c, mpf_t f);
  void mpc_mul_ui (mpc_t rc, mpc_t c, unsigned long int i);
  void mpc_mul_2exp (mpc_t rc, mpc_t c, unsigned long int i);
  void mpc_div (mpc_t rc, mpc_t c1, mpc_t c2);
  void mpc_div_f (mpc_t rc, mpc_t c, mpf_t f);
  void mpc_f_div (mpc_t rc, mpf_t f, mpc_t c);
  void mpc_div_ui (mpc_t rc, mpc_t c, unsigned long int i);
  void mpc_ui_div (mpc_t rc, unsigned long int i, mpc_t c);
  void mpc_div_2exp (mpc_t rc, mpc_t c, unsigned long int i);
  void mpc_pow_si (mpc_t rc, mpc_t c, register signed long int i);
  void mpc_swap (mpc_t c1, mpc_t c2);

/* op= style operators */
  void mpc_smod_eq (mpc_t c);
  void mpc_mod_eq (mpc_t c);
#define mpc_neg_eq(C)           mpc_neg(C, C)
#define mpc_con_eq(C)           mpc_con(C, C)
#define mpc_inv_eq(C)           mpc_inv(C, C)
#define mpc_inv2_eq(C)          mpc_inv2(C, C)
#define mpc_sqr_eq(C)           mpc_sqr(C, C)
  void mpc_rot_eq (mpc_t c);
  void mpc_flip_eq (mpc_t c);
#define mpc_add_eq(R, C)        mpc_add(R, R, C)
#define mpc_add_eq_f(C, R, F)   mpc_add_f(C, C, R, F)
#define mpc_add_eq_ui(C, R, I)  mpc_add_ui(C, C, R, I)
#define mpc_sub_eq(C1, C2)      mpc_sub(C1, C1, C2)
#define mpc_sub_eq_f(C, R, F)   mpc_sub_f(C, C, R, F)
#define mpc_sub_eq_ui(C, R, I)  mpc_sub_ui(C, C, R, I)
#define mpc_ui_sub_eq(C, R, I)  mpc_ui_sub(C, R, I, C)
#define mpc_mul_eq(C1, C2)      mpc_mul(C1, C1, C2)
#define mpc_mul_eq_ui(C, I)     mpc_mul_ui(C, C, I)
#define mpc_mul_eq_f(C, F)      mpc_mul_mpf(C, C, F)
#define mpc_mul_eq_2exp(C, I)   mpc_mul_2exp(C, C, I)
#define mpc_div_eq(C1, C2)      mpc_div(C1, C1, C2)
#define mpc_div2_eq(C1, C2)     mpc_div2(C1, C1, C2)
#define mpc_div_eq_ui(C, I)     mpc_div_ui(C, C, I)
#define mpc_ui_div_eq(C, I)     mpc_ui_div(C, I, C)
#define mpc_div_eq_f(C, F)      mpc_div_f(C, C, F)
#define mpc_div_eq_2exp(C, I)   mpc_div_2exp(C, C, I)
#define mpc_pow_eq_si(C, I)     mpc_pow_si(C, C, I)

/* relational operators */
  int mpc_eq (mpc_t c1, mpc_t c2, unsigned long int i);
  int mpc_eq_zero (mpc_t c);
  int mpc_eq_one (mpc_t c);

/* I/O functions */
  size_t mpc_out_str_2u (FILE * f, int base, size_t n_digits_r,
                         size_t n_digits_i, mpc_t c);
  size_t mpc_out_str_2 (FILE * f, int base, size_t n_digits_r,
                        size_t n_digits_i, mpc_t c);
#define mpc_out_str_u(F, B, D, C)  mpc_out_str_2u(F, B, D, D, C)
#define mpc_out_str(F, B, D, C)  mpc_out_str_2(F, B, D, D, C)
#define mpc_outln_str_u(F, B, D, C)  mpc_out_str_2u(F, B, D, D, C); fputc('\n', F)
#define mpc_outln_str(F, B, D, C)  mpc_out_str_2(F, B, D, D, C); fputc('\n', F)

  size_t mpc_inp_str_u (mpc_t c, FILE * f, int base);
  size_t mpc_inp_str (mpc_t c, FILE * f, int base);

/* vector functions */
#define mpc_valloc(N)           (mpc_t *) malloc((N) * sizeof(mpc_t))
  void mpc_vinit (mpc_t v[], long size);
  void mpc_vinit2 (mpc_t v[], long size, long prec);
  void mpc_vclear (mpc_t v[], long size);
#define mpc_vfree(C)            free(C)

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
