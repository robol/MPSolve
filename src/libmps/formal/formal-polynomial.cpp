#include <mps/mps.h>
#include <exception>
#include <gmpxx.h>
#include <iostream>

using namespace mps::formal;

extern "C" {

  mps_formal_polynomial * 
  mps_formal_polynomial_new_with_monomial (mps_formal_monomial * m)
  {
    Polynomial * p = new Polynomial (*reinterpret_cast<Monomial*> (m));
    return reinterpret_cast<mps_formal_polynomial*> (p);
  }

  mps_formal_polynomial * 
  mps_formal_polynomial_sum_eq (mps_formal_polynomial * p, 
				mps_formal_monomial * m)
  {
    Polynomial * poly = reinterpret_cast<Polynomial*> (p);
    *poly += *reinterpret_cast<Monomial*> (m);
    return reinterpret_cast<mps_formal_polynomial*> (poly);
  }

  mps_formal_polynomial *
  mps_formal_polynomial_sub_eq (mps_formal_polynomial * p, 
				mps_formal_monomial * m)
  {
    Polynomial * poly = reinterpret_cast<Polynomial*> (p);
    *poly -= *reinterpret_cast<Monomial*> (m);
    return reinterpret_cast<mps_formal_polynomial*> (poly);
  }

  mps_formal_polynomial * 
  mps_formal_polynomial_sum_eq_p (mps_formal_polynomial * p, 
				  mps_formal_polynomial * m)
  {
    Polynomial * poly = reinterpret_cast<Polynomial*> (p);
    *poly += *reinterpret_cast<Polynomial*> (m);
    return reinterpret_cast<mps_formal_polynomial*> (poly);
  }

  mps_formal_polynomial *
  mps_formal_polynomial_sub_eq_p (mps_formal_polynomial * p, 
				  mps_formal_polynomial * m)
  {
    Polynomial * poly = reinterpret_cast<Polynomial*> (p);
    *poly -= *reinterpret_cast<Polynomial*> (m);
    return reinterpret_cast<mps_formal_polynomial*> (poly);
  }

  mps_monomial_poly * 
  mps_formal_polynomial_create_monomial_poly (mps_formal_polynomial * p,
					      mps_context * ctx)
  {
    return reinterpret_cast<Polynomial*> (p)->createMonomialPoly (ctx);
  }

  mps_formal_polynomial * 
  mps_formal_polynomial_mul_eq (mps_formal_polynomial * p, 
				mps_formal_polynomial * q)
  {
    Polynomial * poly = reinterpret_cast<Polynomial*> (p);
    *poly *= *reinterpret_cast<Polynomial*>(q);
    return reinterpret_cast<mps_formal_polynomial*> (poly);
  }

  mps_formal_polynomial *
  mps_formal_polynomial_mul (mps_formal_polynomial * p,
			     mps_formal_polynomial * q)
  {
    Polynomial * poly = new Polynomial(*reinterpret_cast<Polynomial*>(p) * 
				       *reinterpret_cast<Polynomial*>(q));
    return reinterpret_cast<mps_formal_polynomial*> (poly);
  }

  void
  mps_formal_polynomial_print (mps_formal_polynomial * p)
  {
    std::cout << *reinterpret_cast<Polynomial*> (p);
  }

  void
  mps_formal_polynomial_free (mps_formal_polynomial * p)
  {
    delete reinterpret_cast<Polynomial*> (p);
  }
  
}

Polynomial::Polynomial()
{
  mMonomials.resize(1);
  mMonomials[0] = Monomial("0", 0);
}

Polynomial::Polynomial(Monomial m)
{
  mMonomials.resize(m.degree() + 1, Monomial("0", 0));
  mMonomials[m.degree()] = m;
}

Polynomial::Polynomial(const Polynomial& rhs)
{
  mMonomials.resize(rhs.degree() + 1, Monomial("0", 0));
  for (int i = 0; i <= rhs.degree(); i++)
    {
      mMonomials[i] = rhs[i];
    }
}

const Monomial
Polynomial::operator[](const int degree) const
{
  if (degree > this->degree() || degree < 0)
    throw std::out_of_range ("Invalid degree specified");
  else
    {
      return mMonomials[degree];
    }
}

Polynomial&
Polynomial::operator-=(const Monomial& m)
{
  Monomial m2 = m;
  *this += (-m2);
  return *this;
}

Polynomial
Polynomial::operator-(const Monomial& m) const
{
  Polynomial out = *this;
  out -= m;
  return out;
}

Polynomial&
Polynomial::operator+=(const Monomial& m)
{ 
  if (m.degree() <= this->degree())
    {
      Monomial currentMonomial = mMonomials[m.degree()];

      if (currentMonomial.isZero())   
	{
	  mMonomials[m.degree()] = m;
	}
      else
	{
	  mMonomials[m.degree()] = Monomial(currentMonomial.coefficientReal() + m.coefficientReal(), 
					    currentMonomial.coefficientImag() + m.coefficientImag(),
					    m.degree());
	}
    }
  else
    {
      mMonomials.resize(m.degree() + 1);
      mMonomials[m.degree()] = m;
    }

  /* Possibly deflate the polynomial, if necessary */
  while (mMonomials[degree()].isZero() && degree() > 0)
    {
      mMonomials.resize(degree());
    }

  return *this;
}

Polynomial
Polynomial::operator+(const Monomial& m) const
{
  Polynomial out = *this;
  out += m;
  return out;
}

Polynomial&
Polynomial::operator+=(const Polynomial& p)
{
  for (int i = 0; i<= p.degree(); i++)
    *this += p[i];
  return *this;
}

Polynomial
Polynomial::operator+(const Polynomial& p) const
{
  Polynomial result = *this;
  result += p;
  return result;
}

Polynomial&
Polynomial::operator-=(const Polynomial& p)
{
  for (int i = 0; i<= p.degree(); i++)
    *this -= p[i];
  return *this;
}

Polynomial
Polynomial::operator-(const Polynomial& p) const
{
  Polynomial result = *this;
  result -= p;
  return result;
}

Polynomial&
Polynomial::operator*=(const Polynomial& other)
{
  Polynomial self = *this * other;

  mMonomials.resize(self.degree() + 1, Monomial("0", 0));
  for (int i = 0; i <= self.degree(); i++)
    {
      mMonomials[i] = self[i];
    }

  return *this;
}

Polynomial
Polynomial::operator*(const Polynomial& other) const
{
  Polynomial result;

  for (int i = 0; i <= degree() + other.degree(); i++)
    {
      for (int j = MAX(0, i - degree()); 
	   j <= MIN(other.degree(), i);
	   j++)
	{
	  result += mMonomials[i-j] * other.mMonomials[j];	  
	}
    }  

  return result;
}

Polynomial::~Polynomial()
{
}

long
Polynomial::degree() const
{
  return mMonomials.size() - 1;
}

Polynomial 
mps::formal::operator+(Monomial a, Monomial b)
{
  Polynomial p(a);
  p += b;
  return p;
}

namespace mps {
  namespace formal {
    std::ostream& 
    operator<<(std::ostream& os, const Polynomial& p)
    {
      os << p[p.degree()];
      
      for (int j = p.degree() - 1; j >= 0; j--)
	{
	  mps::formal::Monomial c = p[j];
	  
	  if (!c.isZero())
	    {
	      if (c.isReal() || c.isImag())
		{
		  if (c.coefficientReal() >= 0 && c.coefficientImag() >= 0)
		    {
		      os << " + " ;
		      os << c;
		    }
		  else
		    {
		      os << " - ";
		      os << -c;
		    }
		}
	      else
		{
		  os << " + " << c;
		}
	    }
	}
      
      return os;
    }
  }
}

mps_monomial_poly *
Polynomial::createMonomialPoly (mps_context * ctx) const
{
  mps_monomial_poly * mp = mps_monomial_poly_new (ctx, degree());
  mpq_t cr, ci;

  mpq_init (cr);
  mpq_init (ci);
  
  for (int i = 0; i <= degree(); i++)
    {
      mpq_set (cr, mMonomials[i].coefficientReal().get_mpq_t());
      mpq_set (ci, mMonomials[i].coefficientImag().get_mpq_t());
      mps_monomial_poly_set_coefficient_q (ctx, mp, i, cr, ci);
    }

  mpq_clear (cr);
  mpq_clear (ci);

  return mp;
}
