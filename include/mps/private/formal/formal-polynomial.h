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

#ifndef MPS_FORMAL_POLYNOMIAL_H_
#define MPS_FORMAL_POLYNOMIAL_H_

MPS_BEGIN_DECLS

MPS_END_DECLS

#ifdef __cplusplus

#include <vector>
#include <stdexcept>
#include <iostream>

namespace mps {
  namespace formal {

    class Polynomial {
    public:
      
      /**
       * @brief Create a default empty Polynomial. 
       */
      Polynomial();

      /**
       * @brief Create a Polynomial made of a single Monomial.
       */
      Polynomial(Monomial m);

      /**
       * @brief Copy constructor for Polynomial. 
       */
      Polynomial(const Polynomial& rhs);

      /**
       * @brief Access a Monomial inside the Polynomial. 
       */
      const Monomial operator[] (const int degree) const;

      /**
       * @brief Add a new Monomial in an existing Polynomial.
       */
      Polynomial& operator+=(const Monomial& m);

      /**
       * @brief Sum two Monomials. 
       */
      Polynomial operator+(const Monomial& m);

      ~Polynomial();

      /**
       * @brief Returns the degree of this polynomial. Empty polynomials
       * have -1 degree. 
       */
      long degree() const;

    private:
      std::vector<Monomial> mMonomials;
    };

    Polynomial operator+(Monomial a, Monomial b);

    std::ostream& operator<<(std::ostream& os, const Polynomial& p);

  }
}

#endif /* ifdef __cplusplus */

#endif /* MPS_FORMAL_POLYNOMIAL_H_ */

