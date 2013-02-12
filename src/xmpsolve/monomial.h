#ifndef XMPSOLVE_MONOMIAL_H
#define XMPSOLVE_MONOMIAL_H

#include <gmp.h>
#include <mps/mps.h>
#include <QString>
#include <QSet>
#include <QChar>

namespace xmpsolve {

/*
 * The parsing in the Monomial class is a rough parser for the following
 * grammar:
 *
 * P -> P | P + M | P - M
 * M -> x | x^E | Cx^E
 * E -> N
 * C -> D | D/D
 * D -> S | Si | (S + Si)
 * S -> G | GeN
 * G -> N | N.N
 */

/**
 * @brief The Monomial class represent a monomial with its coefficient
 * and exponent. It is meant to do the actual parsing work to represent
 * a coefficients given it's string representation.
 *
 * The abstract method addCoefficientToPoly(), in particular, can be
 * used to insert the coefficient into a monomial poly.
 */
class Monomial
{
public:
    /**
     * @param input is the representation of the monomial to parse.
     * @param parent is the parent QObject.
     */
    explicit Monomial(QString input);

    explicit Monomial(double coefficient, int exponent);

    Monomial(const Monomial& other);

    /**
     * @brief isValid checks if the parsing of the input succeded.
     * @return true if the input was valid, and false otherwise.
     */
    bool isValid();

    /**
     * @brief getErrorMessage returns the error in the input string
     * in the case where isValid() is FALSE.
     * @return A QString representation of the error message.
     */
    QString errorMessage();

    /**
     * @brief Degree can be used to access the degree of the monomial
     * @return the degree of the polynomial.
     */
    int degree() const;

    /**
     * @brief addToMonomialPoly can be used to add this this monomial
     * to a given mps_monomial_poly.
     * @param ctx The current mps_context
     * @param poly The poly to which the monomial should be added
     */
    void addToMonomialPoly(mps_context * ctx, mps_monomial_poly * poly);

    /**
     * @brief changeSign Change the sign of the coefficient in front of
     * the monomial.
     */
    void changeSign();

    mpq_t realRationalCoefficient;
    mpq_t imagRationalCoefficient;

    ~Monomial();

    Monomial& operator=(const Monomial& rhs);

private:
    int m_degree;
    QSet<QChar> m_validChars;

    bool m_valid;
    QString m_errorMessage;


    void parseMonomial(QString block);
    void parseCoefficient(QString coefficient);
    void parseNumber(QString number, mpq_t output, mpq_t imag_output);
    void setError(QString message);
    
};

} // namespace xmpsolve

#endif // XMPSOLVE_MONOMIAL_H
