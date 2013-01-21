#include "polynomial.h"

namespace xmpsolve {

Polynomial::Polynomial() :
    m_degree(0)
{
}

// Operator implementation
Polynomial&
Polynomial::operator=(const Polynomial& rhs)
{
    m_degree = rhs.degree();
    return *this;
}

} // namespace xmpsolve
