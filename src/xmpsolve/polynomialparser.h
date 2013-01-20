#ifndef XMPSOLVE_POLYNOMIALPARSER_H
#define XMPSOLVE_POLYNOMIALPARSER_H

#include <QObject>
#include <QList>
#include <mps/mps.h>
#include "monomial.h"

namespace xmpsolve {

class PolynomialParser : public QObject
{
    Q_OBJECT
public:
    explicit PolynomialParser(mps_context * context, QObject *parent = 0);

    /**
     * @brief reset clears all internal data of the parser
     * so it can be re-used.
     */
    void reset();

    /**
     * @brief parse will parse an input polynomial described by the string
     * given as input.
     * @param input is the string representing the input polynomial.
     * @return A pointer to a newly allocated mps_monomial_poly.
     */
    mps_monomial_poly * parse(QString input);

    /**
     * @brief errorMessage can be used to access the last error message in the
     * parser, or the empty string if there was no error.
     * @return A QString with the last error encountered.
     */
    QString errorMessage();

private:
    QList<Monomial*> m_monomials;
    QString m_errorMessage;
    mps_context *m_context;
    int m_degree;
    
signals:
    
public slots:
    
};

} // namespace xmpsolve

#endif // XMPSOLVE_POLYNOMIALPARSER_H
