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

  mps_monomial_poly * 
  mps_formal_polynomial_create_monomial_poly (mps_formal_polynomial * p,
					      mps_context * ctx)
  {
    return reinterpret_cast<Polynomial*> (p)->createMonomialPoly (ctx);
  }

  void
  mps_formal_polynomial_print (mps_formal_polynomial * p)
  {
    std::cout << *reinterpret_cast<Polynomial*> (p);
  }
  
}

Polynomial::Polynomial()
{
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
Polynomial::operator-(const Monomial& m)
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
	  mMonomials[m.degree()] = Monomial(currentMonomial.coefficient() + m.coefficient(), 
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
Polynomial::operator+(const Monomial& m)
{
  Polynomial out = *this;
  out += m;
  return out;
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


std::ostream& 
mps::formal::operator<<(std::ostream& os, const mps::formal::Polynomial& p)
{
  os << p[p.degree()];

  for (int j = p.degree() - 1; j >= 0; j--)
    {
      mps::formal::Monomial c = p[j];

      if (!c.isZero())
	{
	  if (c.coefficient() >= 0)
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
    }

  return os;
}

mps_monomial_poly *
Polynomial::createMonomialPoly (mps_context * ctx) const
{
  mps_monomial_poly * mp = mps_monomial_poly_new (ctx, degree());
  mpq_t zero, c;

  mpq_init (zero);
  mpq_init (c);
  mpq_set_ui (zero, 0U, 1U);
  
  for (int i = 0; i <= degree(); i++)
    {
      mpq_set (c, mMonomials[i].coefficient().get_mpq_t());
      mps_monomial_poly_set_coefficient_q (ctx, mp, i, c, zero);
    }

  mpq_clear (c);
  mpq_clear (zero);

  return mp;
}
