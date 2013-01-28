#ifndef XMPSOLVE_POLYNOMIAL_H
#define XMPSOLVE_POLYNOMIAL_H

#include <QList>
#include "monomial.h"
#include <vector>

namespace xmpsolve {

class Polynomial
{
public:
    explicit Polynomial();
    int degree() const { return m_degree; }

    // Operators
    Polynomial& operator=(const Polynomial& rhs);

    Polynomial& operator+=(const Monomial& rhs);
    Polynomial& operator-=(const Monomial& rhs);
    const Polynomial operator+(const Monomial rhs) const;
    const Polynomial operator-(const Monomial rhs) const;

    /**
     * @brief monomial retrieves the monomial of degree specified
     * @param degree is the degree desired.
     * @return A reference to the monomial of given degree.
     */
    Monomial monomial(int degree) const;

private:
    int m_degree;
    void deflate();

    std::vector<Monomial> m_monomials;
    
};

} // namespace xmpsolve

#endif // XMPSOLVE_POLYNOMIAL_H
