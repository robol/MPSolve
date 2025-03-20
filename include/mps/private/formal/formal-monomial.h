/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
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

struct mps_formal_monomial;

typedef struct mps_formal_monomial mps_formal_monomial;

mps_formal_monomial * mps_formal_monomial_new_with_string (const char *, long);

mps_formal_monomial * mps_formal_monomial_new_with_strings (const char * real, const char * imag, 
							    long degree);

void mps_formal_monomial_free (mps_formal_monomial*);

void mps_formal_monomial_print (mps_formal_monomial*);

mps_formal_monomial * mps_formal_monomial_neg (mps_formal_monomial * m);

mps_formal_monomial * mps_formal_monomial_mul_eq (mps_formal_monomial * m, 
						  mps_formal_monomial * other);

mps_formal_monomial * mps_formal_monomial_mul (mps_formal_monomial * m,
					       mps_formal_monomial * other);

const char * mps_formal_monomial_get_str(mps_formal_monomial * m);

int mps_formal_monomial_degree (mps_formal_monomial *m);

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
       * @brief Create a monomial starting from a string representing the
       * coefficient and an integer representing the degree. 
       */
      Monomial(const char * real_part, const char * imag_part, long degree);

      /**
       * @brief Create a monomial with the given coefficient and of the specified
       * degree. 
       */
      Monomial(const mpq_class coeff, long degree);

      /**
       * @brief Create a monomial with the given coefficient and of the specified
       * degree. 
       */
      Monomial(const mpq_class real, const mpq_class imag, long degree);

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
       * @brief Access the real GMP coefficient stored inside the Monomial.
       */
      const mpq_class coefficientReal() const { return mCoeffR; }
      
      /**
       * @brief Access the imaginary GMP coefficient stored inside the Monomial.
       */
      const mpq_class coefficientImag() const { return mCoeffI; }

      /**
       * @brief Check if this is the zero monomial. 
       */
      bool isZero() const;

      /**
       * @brief Check whether the monomial has a real coefficient. 
       */
      bool isReal() const;

      /**
       * @brief Check whether the monomial has a real coefficient. 
       */
      bool isImag() const;

      Monomial operator-();

      /**
       * @brief Multiply two monomials together. 
       */
      Monomial& operator*=(const Monomial& other);

      /**
       * @brief Multiply two monomials together. 
       */
      Monomial operator*(const Monomial& other) const;

      /**
       * @brief Print a monomial to a stream. 
       */
      friend std::ostream& operator<<(std::ostream& os, const Monomial& l);

    private:
      mpq_class mCoeffR;
      mpq_class mCoeffI;

      long mDegree;
    };

  }
}

#endif /* ifdef __cplusplus */

#endif /* MPS_FORMAL_MONOMIAL_H_ */

