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
  * @brief Implementation in C++ of arithmetic between formal polynomials
  * with rational coefficients. 
  *
  * This implementation is used while parsing polynomials with the lexer 
  * to construct a polynomial incrementally. 
  */

#ifndef MPS_FORMAL_MONOMIAL_H_
#define MPS_FORMAL_MONOMIAL_H_

MPS_BEGIN_DECLS

MPS_END_DECLS

#ifdef __cplusplus

#include <iostream>
#include <gmpxx.h>

namespace mps {
  namespace formal {

    class Monomial {
    public:
      
      /**
       * @brief Create a default Monomial, consisting of the 0 constant with
       * conventional degree -1. 
       */
      Monomial();

      /**
       * @brief Create a monomial starting from a string representing the
       * coefficient and an integer representing the degree. 
       */
      Monomial(const char * coeff_string, long degree);

      /**
       * @brief Create a monomial with the given coefficient and of the specified
       * degree. 
       */
      Monomial(const mpq_class coeff, long degree);

      /**
       * @brief Copy constructor. 
       */
      Monomial(const Monomial& rhs);

      ~Monomial();

      /**
       * @brief Get the degree of the monomial.
       */
      const long degree() const { return mDegree; }

      /**
       * @brief Access the GMP coefficient stored inside the Monomial.
       */
      const mpq_class coefficient() const { return mCoeff; }

      /**
       * @brief Check if this is the zero monomial. 
       */
      bool isZero() const;

      Monomial operator-();

      friend std::ostream& operator<<(std::ostream& os, const Monomial& l);

    private:
      mpq_class mCoeff;
      long mDegree;
    };

    /**
     * @brief A const value with a zeroMonomial that may be used to return const references
     * for methods that need to provide a base case. 
     */
    static const Monomial zeroMonomial("0", 0);

  }
}


#endif /* ifdef __cplusplus */

#endif /* MPS_FORMAL_MONOMIAL_H_ */

