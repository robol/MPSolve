#ifndef XMPSOLVE_MONOMIAL_H
#define XMPSOLVE_MONOMIAL_H

#include <gmp.h>
#include <mps/mps.h>
#include <QString>
#include <QSet>
#include <QChar>

namespace xmpsolve {

enum PolynomialBasis {
    MONOMIAL = 0,
    CHEBYSHEV
};

}

#endif // XMPSOLVE_MONOMIAL_H
