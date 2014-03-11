/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

 /**
  * @file
  * @brief C++ interface to mps_polynomial
  */

#ifndef MPS_POLYNOMIALXX_H_
#define MPS_POLYNOMIALXX_H_

#include <mps/mpsxx.h>

namespace mps {

  class Polynomial : mps_polynomial {

  public:
    /**
     * @brief This constructor has the main role of adjusting the fake vtable in the C
     * struct to reflect the actual content of the C++ implementation that may have been
     * provided in extension to this class. 
     */
    Polynomial ();

    ~Polynomial ();

    /**
     * @brief Public accessor to the degree of the Polynomial. 
     */
    int get_degree () { return degree; }

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
    long int raise_data (mps_context * ctx, long int wp);
    
  public:
    static mps_boolean feval_wrapper (mps_context * ctx, mps_polynomial *p, cplx_t x, cplx_t value, double * error);
    static mps_boolean deval_wrapper (mps_context * ctx, mps_polynomial *p, cdpe_t x, cdpe_t value, rdpe_t error);
    static mps_boolean meval_wrapper (mps_context * ctx, mps_polynomial *p, mpc_t x, mpc_t value, rdpe_t error);
    

  };

}

#endif /* MPS_POLYNOMIALXX_H_ */

