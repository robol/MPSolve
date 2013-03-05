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

#ifndef __MPS_POLYNOMIAL_H
#define __MPS_POLYNOMIAL_H

#include <mps/mps.h>

/* Macro to easy casting of polynomials subclasses */
#define MPS_POLYNOMIAL(t) ((mps_polynomial*) t)

/* A polynomial is completely determined by the functions that allow to operate on it. 
 * The types of these functions is defined in the following. */

/** 
 * @brief The type of a function that evaluates the poynomial (Standard floating point version).
 */
typedef mps_boolean (*mps_polynomial_feval_t) (mps_context  * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error);

/**
 * @brief The type of a function that evaluates the polynomial (CDPE version). 
 */
typedef mps_boolean (*mps_polynomial_deval_t) (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error);

/**
 * @brief The type of a function that evaluates the polynomial (MP version). 
 * The computation must be carried out with the precision of the value \f$x\f$. 
 */
typedef mps_boolean (*mps_polynomial_meval_t) (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error);

/**
 * @brief Function that will be used to deallocate the polynomial on context destruction. 
 */
typedef void (*mps_polynomial_free_t) (mps_context * ctx, mps_polynomial *p);

/**
 * @brief Function that will be used to raise the precision of the coefficients representing
 * the polynomial to the working precision wp. 
 */
typedef long int (*mps_polynomial_raise_data_t) (mps_context * ctx, mps_polynomial * p, long int wp);

/**
 * @brief Function used to determine useful starting approximation.
 *
 * This is the floating point implementation. 
 */
typedef void (*mps_polynomial_fstart_t) (mps_context * ctx, mps_polynomial * p);

/**
 * @brief Function used to determine useful starting approximation.
 *
 * This is the CDPE implementation. 
 */
typedef void (*mps_polynomial_dstart_t) (mps_context * ctx, mps_polynomial * p);

/**
 * @brief Function used to determine useful starting approximation.
 *
 * This is the MP implementation. 
 */
typedef void (*mps_polynomial_mstart_t) (mps_context * ctx, mps_polynomial * p);

/**
 * @brief Function that computes \f$\frac{p}{p'}\f$ (floating point version)
 */
typedef void (*mps_polynomial_fnewton_t) (mps_context * ctx, mps_polynomial * p, 
                                          mps_approximation * root, cplx_t x);
/**
 * @brief Function that computes \f$\frac{p}{p'}\f$ (dpe version)
 */
typedef void (*mps_polynomial_dnewton_t) (mps_context * ctx, mps_polynomial * p, 
                                          mps_approximation * root, cdpe_t corr);

/**
 * @brief Function that computes \f$\frac{p}{p'}\f$ (multiprecision version)
 */
typedef void (*mps_polynomial_mnewton_t) (mps_context * ctx, mps_polynomial * p, 
                                          mps_approximation * root, mpc_t corr);

/** 
 * @brief Function that returns the leading coefficient of the polynomial. 
 * This defaults to the function that returns one (i.e. the default polynomial
 * is monic). 
 */
typedef void (*mps_polynomial_get_leading_coefficient_t) (mps_context * ctx, mps_polynomial * p, mpc_t lc);

/**
 * @brief Struct that represents an abstract polynomial. All the other
 * real polynomial implementations (such as mps_monomial_poly, mps_secular_equation, ...)
 * inherits from this. 
 */
struct mps_polynomial
{
  /**
   * @brief Name of the type. This must be a global static string that
   * can be used to check if a mps_polynomial is of a specific type.
   * It can be NULL to leave the type vague. 
   */
  const char * type_name;

  /**
   * @brief The degree of the polynomial. 
   */
  int degree;

  /**
   * @brief Bits of precision of the coefficients. 
   *
   * The precision used in computation can be adjusted with a call to mps_polynomial_raise_data()
   * but can never be higher than the input precision. 
   */
  long int prec;

  /**
   * @brief Structure of the polynomial, i.e., the algebraic
   * (or non-algebraic) structure where the coefficients
   * are found. 
   */
  mps_structure structure;

  /**
   * @brief Density of the coefficients, or MPS_DENSITY_USER
   * if the coefficients (or the newton fraction) is provided
   * via a user routine
   */
  mps_density density;

  /**
   * @brief This is true if the polynomial has thread-safe methods. Note that
   * this is the default assumption set by mps_polynomial_init(). You should 
   * overwrite after calling it if that's not the case.
   */
  mps_boolean thread_safe;

  /**
   * @brief Method that evaluates the polynomial. 
   */
  mps_polynomial_feval_t feval;

  /**
   * @brief Method that evaluates the polynomial. 
   */
  mps_polynomial_deval_t deval;

  /**
   * @brief Method that evaluates the polynomial. 
   */
  mps_polynomial_meval_t meval;

  /**
   * @brief Method that collocate initial starting points. 
   */
  mps_polynomial_fstart_t fstart;

  /**
   * @brief Method that collocate initial starting points. 
   */
  mps_polynomial_dstart_t dstart;

  /**
   * @brief Method that collocate initial starting points. 
   */
  mps_polynomial_mstart_t mstart;  

  /**
   * @brief Function used to release polynomial resources. 
   */
  mps_polynomial_free_t free;

  /**
   * @brief Function used to raise precision of the coefficients
   * of the representation of the polynomial. 
   */
  mps_polynomial_raise_data_t raise_data;

  /**
   * @brief Function used to compute the Newton correction in a point.
   */
  mps_polynomial_fnewton_t fnewton;

  /**
   * @brief Function used to compute the Newton correction in a point. 
   */
  mps_polynomial_dnewton_t dnewton;

  /**
   * @brief Function used to compute the Newton correction in a point. 
   */
  mps_polynomial_mnewton_t mnewton;

  /**
   * @brief Function used to retrieve the leading coefficient of the
   * polynomial. 
   */
  mps_polynomial_get_leading_coefficient_t get_leading_coefficient;
};

void mps_polynomial_init (mps_context * ctx, mps_polynomial * p);

mps_polynomial * mps_polynomial_new (mps_context * ctx);

mps_boolean mps_polynomial_check_type (mps_polynomial * p, const char * type_name);

mps_boolean mps_polynomial_feval (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error);

mps_boolean mps_polynomial_deval (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error);

mps_boolean mps_polynomial_meval (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error);

void mps_polynomial_fstart (mps_context * ctx, mps_polynomial * p);

void mps_polynomial_dstart (mps_context * ctx, mps_polynomial * p);

void mps_polynomial_mstart (mps_context * ctx, mps_polynomial * p);

void mps_polynomial_free (mps_context * ctx, mps_polynomial * p);

void mps_polynomial_fnewton (mps_context * ctx, mps_polynomial *p, 
                             mps_approximation * root, cplx_t corr);

void mps_polynomial_dnewton (mps_context * ctx, mps_polynomial *p, 
                             mps_approximation * root, cdpe_t corr);

void mps_polynomial_mnewton (mps_context * ctx, mps_polynomial *p, 
                             mps_approximation * root, mpc_t corr);

void mps_polynomial_get_leading_coefficient (mps_context * ctx, mps_polynomial * p, mpc_t lc);

long int mps_polynomial_raise_data (mps_context * ctx, mps_polynomial * p, long int wp);

#ifdef _MPS_PRIVATE
  /* Private implementation details */
#endif

#endif
