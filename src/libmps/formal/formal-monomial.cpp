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
  
  mps_formal_monomial *
  mps_formal_monomial_new_with_strings (const char * real, const char * imag, 
					long degree)
  {
    Monomial * m = new Monomial (real, imag, degree);
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

  mps_formal_monomial * 
  mps_formal_monomial_mul_eq (mps_formal_monomial * m, 
			      mps_formal_monomial * other)
  {
    Monomial * mm = reinterpret_cast<Monomial*> (m);
    *mm *= *reinterpret_cast<Monomial*> (other);
    return reinterpret_cast<mps_formal_monomial*> (mm);
  }

  mps_formal_monomial * 
  mps_formal_monomial_mul (mps_formal_monomial * m,
			   mps_formal_monomial * other)
  {
    Monomial * result = new Monomial (*reinterpret_cast<Monomial*> (m) * 
				      *reinterpret_cast<Monomial*> (other));
    return reinterpret_cast<mps_formal_monomial*> (result);
  }

}

Monomial::Monomial()
{
  mCoeffR = 0;
  mCoeffI = 0;
  mDegree = 0;
}

Monomial::Monomial(const char * coeff_string, long degree)
{
  char * er = mps_utils_build_equivalent_rational_string (NULL, coeff_string);

  mCoeffR = er;
  mDegree = degree;
  free (er);
}

Monomial::Monomial(const char * real_part, const char * imag_part, long degree)
{
  char * er = mps_utils_build_equivalent_rational_string (NULL, real_part);
  char * ei = mps_utils_build_equivalent_rational_string (NULL, imag_part);

  mDegree = degree;

  mCoeffR = er;
  mCoeffI = ei;

  free (er);
  free (ei);
}

Monomial::Monomial(const mpq_class coeff, long degree)
{
  mCoeffR = coeff;
  mCoeffR.canonicalize();
  mDegree = degree;
}

Monomial::Monomial(const mpq_class realpart, const mpq_class imagpart, long degree)
{
  mCoeffR = realpart;
  mCoeffI = imagpart;

  mCoeffR.canonicalize();
  mCoeffI.canonicalize();

  mDegree = degree;
}

Monomial::Monomial(const Monomial& rhs)
{
  mCoeffR = rhs.coefficientReal();
  mCoeffI = rhs.coefficientImag();
  mDegree = rhs.degree();
}

Monomial::~Monomial()
{
}

bool
Monomial::isZero() const
{
  return mCoeffR == 0 && mCoeffI == 0;
}

bool
Monomial::isReal() const
{
  return mCoeffI == 0;
}

bool
Monomial::isImag() const
{
  return mCoeffR == 0;
}



Monomial
Monomial::operator-()
{
  return Monomial(-mCoeffR, -mCoeffI, mDegree);
}

Monomial& 
Monomial::operator*=(const Monomial& other)
{
  mpq_class tmp;

  tmp = mCoeffR * other.mCoeffR - mCoeffI * other.mCoeffI;
  mCoeffI = mCoeffI * other.mCoeffR + mCoeffR * other.mCoeffI;
  mCoeffR = tmp;
  mDegree += other.mDegree;

  return *this;
}

Monomial 
Monomial::operator*(const Monomial& other) const
{
  Monomial result = *this;
  result *= other;
  return result;
}


std::ostream&
mps::formal::operator<<(std::ostream& os, const mps::formal::Monomial& m)
{
  if (m.isReal())
    os << m.mCoeffR;
  else if (m.mCoeffR == 0)
    os << m.mCoeffI << "i";
  else 
    os << "(" << m.mCoeffR << (m.mCoeffI > 0 ? "+" : "-")
       << (m.mCoeffI > 0 ? m.mCoeffI : -m.mCoeffI) << "i)";

  switch (m.degree())
    {
    case 0:
      break;
    case 1:
      os << "x";
      break;
    default:
      os << "x^" << m.mDegree;
      break;
    }

  return os;
}

