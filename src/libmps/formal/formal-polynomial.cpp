#include <mps/mps.h>
#include <exception>
#include <gmpxx.h>

#include <iostream>

using namespace mps::formal;

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

  return os;
}
