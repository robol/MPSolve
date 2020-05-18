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

/*! @cond PRIVATE */
typedef struct {
  mpf_t r, i;
} __mpcf_struct;
/*! @endcond */

typedef __mpcf_struct mpcf_t[1];

/***********************************************************
**            macros for fields access                    **
***********************************************************/

/*
 * #define mpcf_Val(C)             (*(C))
 */
#define mpcf_Re(C)               ((C)->r)
#define mpcf_Im(C)               ((C)->i)
#define mpcf_Addr(C)             ((__mpcf_struct*)C)
#define mpcf_Move(C1, C2)        (*(C1) = *(C2))

/***********************************************************
**            mpcf_t functions                             **
***********************************************************/

/* constructors/destructors */
void mpcf_init (mpcf_t c);
void mpcf_init2 (mpcf_t c, unsigned long int prec);
void mpcf_clear (mpcf_t c);

void mpcf_set_prec (mpcf_t c, unsigned long int prec);
unsigned long int mpcf_get_prec (const mpcf_t c);
void mpcf_set_prec_raw (mpcf_t c, unsigned long int prec);

/* initializers */
void mpcf_set (mpcf_t rc, const mpcf_t c);
void mpcf_set_ui (mpcf_t c, unsigned long int ir, unsigned long int ii);
void mpcf_set_si (mpcf_t c, signed long int ir, signed long int ii);
void mpcf_set_d (mpcf_t c, double dr, double di);
void mpcf_set_z (mpcf_t c, mpz_t zr, mpz_t zi);
void mpcf_set_q (mpcf_t c, mpq_t qr, mpq_t qi);
void mpcf_set_f (mpcf_t c, mpf_t fr, mpf_t fi);
int mpcf_set_str (mpcf_t c, char *sr, char *si, int base);

void mpcf_init_set (mpcf_t rc, mpcf_t c);
void mpcf_init_set_ui (mpcf_t c, unsigned long int ir, unsigned long int ii);
void mpcf_init_set_si (mpcf_t c, signed long int ir, signed long int ii);
void mpcf_init_set_d (mpcf_t c, double dr, double di);
void mpcf_init_set_f (mpcf_t c, mpf_t fr, mpf_t fi);
int mpcf_init_set_str (mpcf_t c, char *sr, char *si, int base);

/* unary functions */
void mpcf_neg (mpcf_t rc, mpcf_t c);
void mpcf_smod (mpf_t f, mpcf_t c);
void mpcf_rmod (rdpe_t r, mpcf_t c);
void mpcf_mod (mpf_t f, mpcf_t c);
void mpcf_con (mpcf_t rc, mpcf_t c);
void mpcf_inv (mpcf_t rc, mpcf_t c);
void mpcf_inv2 (mpcf_t rc, mpcf_t c);
void mpcf_sqr (mpcf_t rc, mpcf_t c);
void mpcf_rot (mpcf_t rc, mpcf_t c);
void mpcf_flip (mpcf_t rc, mpcf_t c);

/* binary functions */
void mpcf_add (mpcf_t rc, mpcf_t c1, mpcf_t c2);
void mpcf_add_f (mpcf_t rc, mpcf_t c, mpf_t f);
void mpcf_add_ui (mpcf_t rc, mpcf_t c, unsigned long int r,
                 unsigned long int i);
void mpcf_sub (mpcf_t rc, mpcf_t c1, mpcf_t c2);
void mpcf_sub_f (mpcf_t rc, mpcf_t c, mpf_t f);
void mpcf_f_sub (mpcf_t rc, mpf_t f, mpcf_t c);
void mpcf_sub_ui (mpcf_t rc, mpcf_t c, unsigned long int r,
                 unsigned long int i);
void mpcf_ui_sub (mpcf_t rc, unsigned long int r, unsigned long int i,
                 mpcf_t c);
void mpcf_mul (mpcf_t rc, mpcf_t c1, mpcf_t c2);
void mpcf_mul_f (mpcf_t rc, mpcf_t c, mpf_t f);
void mpcf_mul_ui (mpcf_t rc, mpcf_t c, unsigned long int i);
void mpcf_mul_2exp (mpcf_t rc, mpcf_t c, unsigned long int i);
void mpcf_div (mpcf_t rc, mpcf_t c1, mpcf_t c2);
void mpcf_div_f (mpcf_t rc, mpcf_t c, mpf_t f);
void mpcf_f_div (mpcf_t rc, mpf_t f, mpcf_t c);
void mpcf_div_ui (mpcf_t rc, mpcf_t c, unsigned long int i);
void mpcf_ui_div (mpcf_t rc, unsigned long int i, mpcf_t c);
void mpcf_div_2exp (mpcf_t rc, mpcf_t c, unsigned long int i);
void mpcf_pow_si (mpcf_t rc, mpcf_t c, register signed long int i);
void mpcf_swap (mpcf_t c1, mpcf_t c2);

/* op= style operators */
void mpcf_smod_eq (mpcf_t c);
void mpcf_mod_eq (mpcf_t c);
#define mpcf_neg_eq(C)           mpcf_neg (C, C)
#define mpcf_con_eq(C)           mpcf_con (C, C)
#define mpcf_inv_eq(C)           mpcf_inv (C, C)
#define mpcf_inv2_eq(C)          mpcf_inv2 (C, C)
#define mpcf_sqr_eq(C)           mpcf_sqr (C, C)
void mpcf_rot_eq (mpcf_t c);
void mpcf_flip_eq (mpcf_t c);
#define mpcf_add_eq(R, C)        mpcf_add (R, R, C)
#define mpcf_add_eq_f(C, F)      mpcf_add_f (C, C, F)
#define mpcf_add_eq_ui(C, R, I)  mpcf_add_ui (C, C, R, I)
#define mpcf_sub_eq(C1, C2)      mpcf_sub (C1, C1, C2)
#define mpcf_sub_eq_f(C, R, F)   mpcf_sub_f (C, C, R, F)
#define mpcf_sub_eq_ui(C, R, I)  mpcf_sub_ui (C, C, R, I)
#define mpcf_ui_sub_eq(C, R, I)  mpcf_ui_sub (C, R, I, C)
#define mpcf_mul_eq(C1, C2)      mpcf_mul (C1, C1, C2)
#define mpcf_mul_eq_ui(C, I)     mpcf_mul_ui (C, C, I)
#define mpcf_mul_eq_f(C, F)      mpcf_mul_mpf (C, C, F)
#define mpcf_mul_eq_2exp(C, I)   mpcf_mul_2exp (C, C, I)
#define mpcf_div_eq(C1, C2)      mpcf_div (C1, C1, C2)
#define mpcf_div2_eq(C1, C2)     mpcf_div2 (C1, C1, C2)
#define mpcf_div_eq_ui(C, I)     mpcf_div_ui (C, C, I)
#define mpcf_ui_div_eq(C, I)     mpcf_ui_div (C, I, C)
#define mpcf_div_eq_f(C, F)      mpcf_div_f (C, C, F)
#define mpcf_div_eq_2exp(C, I)   mpcf_div_2exp (C, C, I)
#define mpcf_pow_eq_si(C, I)     mpcf_pow_si (C, C, I)

/* relational operators */
int mpcf_eq (mpcf_t c1, mpcf_t c2, unsigned long int i);
int mpcf_eq_zero (mpcf_t c);
int mpcf_eq_one (mpcf_t c);

/* I/O functions */
size_t mpcf_out_str_2u (FILE * f, int base, size_t n_digits_r,
                       size_t n_digits_i, mpcf_t c);
size_t mpcf_out_str_2 (FILE * f, int base, size_t n_digits_r,
                      size_t n_digits_i, mpcf_t c);
#define mpcf_out_str_u(F, B, D, C)  mpcf_out_str_2u (F, B, D, D, C)
#define mpcf_out_str(F, B, D, C)  mpcf_out_str_2 (F, B, D, D, C)
#define mpcf_outln_str_u(F, B, D, C)  mpcf_out_str_2u (F, B, D, D, C); fputc ('\n', F)
#define mpcf_outln_str(F, B, D, C)  mpcf_out_str_2 (F, B, D, D, C); fputc ('\n', F)

size_t mpcf_inp_str_u (mpcf_t c, FILE * f, int base);
size_t mpcf_inp_str (mpcf_t c, FILE * f, int base);

/* vector functions */
#define mpcf_valloc(N)           (mpcf_t*)malloc ((N)*sizeof(mpcf_t))
void mpcf_vinit (mpcf_t v[], long size);
void mpcf_vinit2 (mpcf_t v[], long size, long prec);
void mpcf_vclear (mpcf_t v[], long size);
#define mpcf_vfree(C)            free (C)

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
