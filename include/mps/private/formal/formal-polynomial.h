/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
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

struct mps_formal_polynomial;

typedef struct mps_formal_polynomial mps_formal_polynomial;

mps_formal_polynomial * mps_formal_polynomial_new_with_monomial (mps_formal_monomial *);

mps_formal_polynomial * mps_formal_polynomial_sum_eq (mps_formal_polynomial * p, 
						      mps_formal_monomial * m);

mps_formal_polynomial * mps_formal_polynomial_sub_eq (mps_formal_polynomial * p, 
						      mps_formal_monomial * m);

mps_formal_polynomial * mps_formal_polynomial_sum_eq_p (mps_formal_polynomial * p, 
							mps_formal_polynomial * m);

mps_formal_polynomial * mps_formal_polynomial_sub_eq_p (mps_formal_polynomial * p, 
							mps_formal_polynomial * m);

mps_monomial_poly * mps_formal_polynomial_create_monomial_poly (mps_formal_polynomial * p,
								mps_context * ctx);

mps_formal_polynomial * mps_formal_polynomial_mul (mps_formal_polynomial * p,
						   mps_formal_polynomial * q);

mps_formal_polynomial * mps_formal_polynomial_mul_eq (mps_formal_polynomial * p,
						      mps_formal_polynomial * q);

void mps_formal_polynomial_print (mps_formal_polynomial * p);

void mps_formal_polynomial_free (mps_formal_polynomial * p);

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
      Polynomial operator+(const Monomial& m) const;

      /**
       * @brief Add a new Polynomial in an existing Polynomial.
       */
      Polynomial& operator+=(const Polynomial& p);

      /**
       * @brief Sum two Polynomials. 
       */
      Polynomial operator+(const Polynomial& p) const;

      /**
       * @brief Subtract a new Polynomial from an existing Polynomial.
       */
      Polynomial& operator-=(const Polynomial& p);

      /**
       * @brief Subtract two Polynomials. 
       */
      Polynomial operator-(const Polynomial& p) const;

      /**
       * @brief Subtract a new Monomial in an existing Polynomial.
       */
      Polynomial& operator-=(const Monomial& m);

      /**
       * @brief Subtract two Monomials. 
       */
      Polynomial operator-(const Monomial& m) const;

      /**
       * @brief Multiply two polynomials together. 
       */
      Polynomial& operator*=(const Polynomial& other);

      /**
       * @brief Multiply two polynomials together. 
       */
      Polynomial operator*(const Polynomial& other) const;

      ~Polynomial();

      /**
       * @brief Returns the degree of this polynomial. Empty polynomials
       * have -1 degree. 
       */
      long degree() const;

      /**
       * @brief Create a representation of this polynomial as a mps_monomial_poly. 
       */
      mps_monomial_poly * createMonomialPoly (mps_context * ctx) const;

    private:
      std::vector<Monomial> mMonomials;
    };

    Polynomial operator+(Monomial a, Monomial b);

    std::ostream& operator<<(std::ostream& os, const Polynomial& p);

  }
}

#endif /* ifdef __cplusplus */

#endif /* MPS_FORMAL_POLYNOMIAL_H_ */

