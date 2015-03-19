/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#ifndef MPS_POLYNOMIAL_H_
#define MPS_POLYNOMIAL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <mps/mps.h>

/* Macro that can be used to enforce a sort of type-safe casting between
 * mps_polynomial "subclasses". Please note that this does not guarantee
 * type-safeness at all if you cast other types. */
#define MPS_POLYNOMIAL_CAST(typename, t) ((typename*)(mps_polynomial_cast (# typename, (mps_polynomial*)t)))

/* Macro to ease casting of polynomials subclasses */
#define MPS_POLYNOMIAL(t) (MPS_POLYNOMIAL_CAST (mps_polynomial, t))

/* A polynomial is completely determined by the functions that allow to operate on it.
 * The types of these functions is defined in the following. */

/**
 * @brief The type of a function that evaluates the poynomial (Standard floating point version).
 */
typedef mps_boolean (*mps_polynomial_feval_t)(mps_context  * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error);

/**
 * @brief The type of a function that evaluates the polynomial (CDPE version).
 */
typedef mps_boolean (*mps_polynomial_deval_t)(mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error);

/**
 * @brief The type of a function that evaluates the polynomial (MP version).
 * The computation must be carried out with the precision of the value \f$x\f$.
 */
typedef mps_boolean (*mps_polynomial_meval_t)(mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error);

/**
 * @brief Function that will be used to deallocate the polynomial on context destruction.
 */
typedef void (*mps_polynomial_free_t)(mps_context * ctx, mps_polynomial *p);

/**
 * @brief Function that will be used to raise the precision of the coefficients representing
 * the polynomial to the working precision wp.
 */
typedef long int (*mps_polynomial_raise_data_t)(mps_context * ctx, mps_polynomial * p, long int wp);

/**
 * @brief Function used to determine useful starting approximation.
 *
 * This is the floating point implementation.
 */
typedef void (*mps_polynomial_fstart_t)(mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);

/**
 * @brief Function used to determine useful starting approximation.
 *
 * This is the CDPE implementation.
 */
typedef void (*mps_polynomial_dstart_t)(mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);

/**
 * @brief Function used to determine useful starting approximation.
 *
 * This is the MP implementation.
 */
typedef void (*mps_polynomial_mstart_t)(mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);

/**
 * @brief Function that computes \f$\frac{p}{p'}\f$ (floating point version)
 */
typedef void (*mps_polynomial_fnewton_t)(mps_context * ctx, mps_polynomial * p,
                                         mps_approximation * root, cplx_t x);
/**
 * @brief Function that computes \f$\frac{p}{p'}\f$ (dpe version)
 */
typedef void (*mps_polynomial_dnewton_t)(mps_context * ctx, mps_polynomial * p,
                                         mps_approximation * root, cdpe_t corr);

/**
 * @brief Function that computes \f$\frac{p}{p'}\f$ (multiprecision version)
 */
typedef void (*mps_polynomial_mnewton_t)(mps_context * ctx, mps_polynomial * p,
                                         mps_approximation * root, mpc_t corr, long int wp);

/**
 * @brief Function that returns the leading coefficient of the polynomial.
 * This defaults to the function that returns one (i.e. the default polynomial
 * is monic).
 */
typedef void (*mps_polynomial_get_leading_coefficient_t)(mps_context * ctx, mps_polynomial * p, mpc_t lc);

/**
 * @brief Struct that represents an abstract polynomial. All the other
 * real polynomial implementations (such as mps_monomial_poly, mps_secular_equation, ...)
 * inherits from this.
 */
struct mps_polynomial {
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

mps_polynomial * mps_polynomial_cast (const char *type_name, mps_polynomial *p);

mps_boolean mps_polynomial_feval (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error);

mps_boolean mps_polynomial_deval (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error);

mps_boolean mps_polynomial_meval (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error);

void mps_polynomial_fstart (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);

void mps_polynomial_dstart (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);

void mps_polynomial_mstart (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);

void mps_polynomial_free (mps_context * ctx, mps_polynomial * p);

void mps_polynomial_fnewton (mps_context * ctx, mps_polynomial *p,
                             mps_approximation * root, cplx_t corr);

void mps_polynomial_dnewton (mps_context * ctx, mps_polynomial *p,
                             mps_approximation * root, cdpe_t corr);

void mps_polynomial_mnewton (mps_context * ctx, mps_polynomial *p,
                             mps_approximation * root, mpc_t corr, long int wp);

void mps_polynomial_get_leading_coefficient (mps_context * ctx, mps_polynomial * p, mpc_t lc);

long int mps_polynomial_raise_data (mps_context * ctx, mps_polynomial * p, long int wp);

void mps_polynomial_set_input_prec (mps_context * ctx, mps_polynomial * p, long int prec);

/* functions in general-starting.c */
void mps_general_fstart (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);
void mps_general_dstart (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);
void mps_general_mstart (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);

#ifdef _MPS_PRIVATE
/* Private implementation details */
#endif

#ifdef __cplusplus
}
#endif


/*
 * C++ wrapper around mps_polynomial.
 */

#ifdef __cplusplus

namespace mps {
  class Polynomial : public mps_polynomial {
public:
    /**
     * @brief This constructor has the main role of adjusting the fake vtable in the C
     * struct to reflect the actual content of the C++ implementation that may have been
     * provided in extension to this class.
     */
    explicit Polynomial (mps_context * ctx, const char * type_name = "mps_polynomial");

    virtual ~Polynomial ();

    /**
     * @brief Public accessor to the degree of the Polynomial.
     */
    int get_degree ()
    {
      return degree;
    }

    /**
     * @brief Evaluate the polynomial at a point.
     *
     * This method should be overloaded by subclasses of {@link Polynomial} in order
     * to provide the necessary methods to MPSolve.
     *
     * @param x The point where the {@link Polynomial} should be evaluted.
     * @param value The storage where the result of the evaluation will be stored.
     * @param error An upper bound to the error that has been computed in this operation.
     *
     * @return true if the operation was successful, false in case an exception has been
     * encountered.
     */
    virtual bool eval (mps_context * ctx, cplx_t x, cplx_t value, double * error) = 0;

    /**
     * @brief Evaluate the polynomial at a point.
     *
     * This method should be overloaded by subclasses of {@link Polynomial} in order
     * to provide the necessary methods to MPSolve.
     *
     * @param x The point where the {@link Polynomial} should be evaluted.
     * @param value The storage where the result of the evaluation will be stored.
     * @param error An upper bound to the error that has been computed in this operation.
     *
     * @return true if the operation was successful, false in case an exception has been
     * encountered.
     */
    virtual bool eval (mps_context * ctx, cdpe_t x, cdpe_t value, rdpe_t error) = 0;

    /**
     * @brief Evaluate the polynomial at a point.
     *
     * This method should be overloaded by subclasses of {@link Polynomial} in order
     * to provide the necessary methods to MPSolve.
     *
     * @param x The point where the {@link Polynomial} should be evaluted.
     * @param value The storage where the result of the evaluation will be stored.
     * @param error An upper bound to the error that has been computed in this operation.
     *
     * @return true if the operation was successful, false in case an exception has been
     * encountered.
     */
    virtual bool eval (mps_context * ctx, mpc_t x, mpc_t value, rdpe_t error) = 0;

    /**
     * @brief Raise the working precision of this polynomial to the specified
     * value.
     *
     * Note that this might be a no-op on polynomials that are defined implicitly or
     * without the need for explicit coefficients.
     */
    virtual long int raise_data_wp (mps_context * ctx, long int wp);

    virtual void start_fp (mps_context * ctx, mps_approximation ** approximations);
    virtual void start_dpe (mps_context * ctx, mps_approximation ** approximations);
    virtual void start_mp (mps_context * ctx, mps_approximation ** approximations);

    virtual void get_leading_coefficient (mps_context * ctx, mpc_t lc);

    virtual void newton (mps_context * ctx, mps_approximation * a, cplx_t x) = 0;
    virtual void newton (mps_context * ctx, mps_approximation * a, cdpe_t x) = 0;
    virtual void newton (mps_context * ctx, mps_approximation * a, mpc_t x, long int wp) = 0;

public:
    static mps_boolean feval_wrapper (mps_context * ctx, mps_polynomial *p,
                                      cplx_t x, cplx_t value, double * error);

    static mps_boolean deval_wrapper (mps_context * ctx, mps_polynomial *p,
                                      cdpe_t x, cdpe_t value, rdpe_t error);

    static mps_boolean meval_wrapper (mps_context * ctx, mps_polynomial *p,
                                      mpc_t x, mpc_t value, rdpe_t error);

    static void free_wrapper (mps_context * ctx, mps_polynomial * p);

    static long int raise_data_wrapper (mps_context * ctx, mps_polynomial * p,
                                        long int wp);

    static void fstart_wrapper (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);
    static void dstart_wrapper (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);
    static void mstart_wrapper (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations);

    static void fnewton_wrapper (mps_context * ctx, mps_polynomial * p,
                                 mps_approximation * a, cplx_t x);
    static void dnewton_wrapper (mps_context * ctx, mps_polynomial * p,
                                 mps_approximation * a, cdpe_t x);
    static void mnewton_wrapper (mps_context * ctx, mps_polynomial * p,
                                 mps_approximation * a, mpc_t x,
                                 long int wp);

    static void get_leading_coefficient_wrapper (mps_context * ctx, mps_polynomial * p,
                                                 mpc_t leading_coefficient);
  };
}

#endif /* __cplusplus */

#endif
