#include "polynomial.h"

namespace xmpsolve {

Polynomial::Polynomial() :
    m_degree(0)
{
    // We are the zero polynomial by default.
    m_monomials.push_back(Monomial("0"));
}

// Operator implementation
Polynomial&
Polynomial::operator=(const Polynomial& rhs)
{
    m_degree = rhs.degree();

    for(int i = 0; i <= m_degree; i++)
    {
        m_monomials.push_back(rhs.monomial(i));
    }

    return *this;
}

Monomial
Polynomial::monomial(int degree) const
{
    if (degree > m_degree)
        return Monomial(0, degree);

    return m_monomials[degree];
}

void
Polynomial::deflate()
{
    mpq_t zero;
    mpq_init(zero);
    mpq_set_ui(zero, 0U, 0U);
    bool deflatable = true;

    for(int i = m_degree; i > 0 && deflatable; i--)
    {
        Monomial m = monomial(i);
        if (mpq_equal (m.realRationalCoefficient, zero) &&
                mpq_equal(m.imagRationalCoefficient, zero)) {
            m_degree = m_degree - 1;
            m_monomials.pop_back();
        }
        else {
            deflatable = false;
        }
    }

    mpq_clear(zero);
}

Polynomial&
Polynomial::operator+=(const Monomial& rhs)
{
    Monomial lhs = monomial(rhs.degree());

    mpq_add(lhs.realRationalCoefficient, lhs.realRationalCoefficient, rhs.realRationalCoefficient);
    mpq_add(lhs.imagRationalCoefficient, lhs.imagRationalCoefficient, rhs.imagRationalCoefficient);

    // Fill with zero monomials the empty fields
    if (rhs.degree() > m_degree) {
        for (int i = m_degree + 1; i <= rhs.degree(); i++) {
            m_monomials.push_back(Monomial(0, i));
        }
        m_degree = rhs.degree();
    }

    m_monomials.at(rhs.degree()) = lhs;
    deflate();

    return *this;
}

Polynomial&
Polynomial::operator-=(const Monomial& rhs)
{
    Monomial lhs = monomial(rhs.degree());

    mpq_sub(lhs.realRationalCoefficient, lhs.realRationalCoefficient, rhs.realRationalCoefficient);
    mpq_sub(lhs.imagRationalCoefficient, lhs.imagRationalCoefficient, rhs.imagRationalCoefficient);

    // Fill with zero monomials the empty fields
    if (rhs.degree() > m_degree) {
        for (int i = m_degree + 1; i < rhs.degree(); i++) {
            m_monomials.push_back(Monomial(0, i));
        }
        m_degree = rhs.degree();
    }

    m_monomials.at(rhs.degree()) = lhs;
    deflate();

    return *this;
}

const Polynomial
Polynomial::operator+(const Monomial rhs) const
{
    Polynomial result = *this;
    result += rhs;
    return result;
}

const Polynomial
Polynomial::operator-(const Monomial rhs) const
{
    Polynomial result = *this;
    result -= rhs;
    return result;
}

} // namespace xmpsolve
