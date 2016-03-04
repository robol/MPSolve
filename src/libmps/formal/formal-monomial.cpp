#include <mps/mps.h>
#include <iostream>

using namespace mps::formal;

extern "C" {

mps_formal_monomial * 
mps_formal_monomial_new_with_string (const char * coeff_string, long degree)
{
  Monomial * m = new Monomial (coeff_string, degree);
  return reinterpret_cast<mps_formal_monomial*> (m);
}

void
mps_formal_monomial_free (mps_formal_monomial* m)
{
  delete reinterpret_cast<Monomial*> (m);
}

  void
  mps_formal_monomial_print (mps_formal_monomial * m)
  {
    std::cout << *reinterpret_cast<Monomial*> (m);
  }

  mps_formal_monomial *
  mps_formal_monomial_neg (mps_formal_monomial * m)
  {
    Monomial *m2 = new Monomial (- *reinterpret_cast<Monomial*> (m));
    return reinterpret_cast<mps_formal_monomial*> (m2);
  }

}

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
  mCoeff.canonicalize();
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

