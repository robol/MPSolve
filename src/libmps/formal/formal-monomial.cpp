#include <mps/mps.h>
#include <iostream>

using namespace mps::formal;

Monomial::Monomial()
{
  mCoeff = 1;
  mDegree = -1;
}

Monomial::Monomial(const char * coeff_string, long degree)
{
  char * er = mps_utils_build_equivalent_rational_string (NULL, coeff_string);
  mCoeff = er;
  mDegree = degree;
  free (er);
}


Monomial::Monomial(const mpq_class coeff, long degree)
{
  mCoeff = coeff;
  mDegree = degree;
}

Monomial::Monomial(const Monomial& rhs)
{
  mCoeff = rhs.coefficient();
  mDegree = rhs.degree();
}

Monomial::~Monomial()
{
}

bool
Monomial::isZero() const
{
  return mCoeff == 0;
}

Monomial
Monomial::operator-()
{
  return Monomial(-mCoeff, mDegree);
}

std::ostream&
mps::formal::operator<<(std::ostream& os, const mps::formal::Monomial& m)
{
  switch (m.degree())
    {
    case 0:
      os << m.mCoeff;
      break;
    case 1:
      os << m.mCoeff << "x";
      break;
    default:
      os << m.mCoeff << "x^" << m.mDegree;
      break;
    }

  return os;
}
