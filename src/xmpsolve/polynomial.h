#ifndef XMPSOLVE_POLYNOMIAL_H
#define XMPSOLVE_POLYNOMIAL_H

#include <QList>
#include "monomial.h"

namespace xmpsolve {

class Polynomial
{
public:
    explicit Polynomial();
    int degree() const { return m_degree; }

    // Operators
    Polynomial& operator=(const Polynomial& rhs);


private:
    int m_degree;
    QList<Monomial> m_monomials;
    
};

} // namespace xmpsolve

#endif // XMPSOLVE_POLYNOMIAL_H
