#include "monomial.h"
#include <QDebug>
#include <QSet>

namespace xmpsolve {

Monomial::Monomial(QString input)
{
    mpq_init(realRationalCoefficient);
    mpq_init(imagRationalCoefficient);

    m_valid = true;
    m_validChars << '1' << '2' << '3' << '4' << '5' << '6' << '7' << '8'
                 << '9' << '0' << 'i' << 'I' << '(' << ')' << '/' << '+'
                 << '-' << 'e' << '.';

    // Start the parsing
    parseMonomial(input.trimmed());

    // Canonicalize the rational number
    mpq_canonicalize(realRationalCoefficient);
    mpq_canonicalize(imagRationalCoefficient);
}

Monomial::Monomial(double coefficient, int exponent)
{
    mpq_init(realRationalCoefficient);
    mpq_init(imagRationalCoefficient);

    mpq_set_d(realRationalCoefficient, coefficient);
    mpq_set_si(imagRationalCoefficient, 0, 0);

    m_degree = exponent;

    // Canonicalize the rational number
    mpq_canonicalize(realRationalCoefficient);
    mpq_canonicalize(imagRationalCoefficient);
}

Monomial::Monomial(const Monomial &other)
{
    mpq_init(realRationalCoefficient);
    mpq_init(imagRationalCoefficient);

    *this = other;
}

void
Monomial::parseMonomial(QString input)
{
    bool conversionOk;

    if (input == QString("x")) {
        m_degree = 1;
        mpq_set_str (realRationalCoefficient, "1", 10);
        mpq_set_str (imagRationalCoefficient, "0", 10);
        return;
    }

    if (input.startsWith("x")) {
        // Check that we really have exponentiation here.
        if (input.at(1) != '^') {
            setError(QObject::QObject::tr("unexpected character %1, expected ^").arg(input.at(1)));
            return;
        }

        // Check which is the degree of this term
        QString exponent = input.right(input.length() - 2);
        m_degree = exponent.toInt(&conversionOk);

        if (!conversionOk) {
            setError(QObject::QObject::tr("not valid exponent: %1").arg(exponent));
            return;
        }

        mpq_set_str(realRationalCoefficient, "1", 10);
        mpq_set_str(imagRationalCoefficient, "0", 10);
    }

    // In the last case we have a monomial of the type Cx^E
    int exp_pos = input.indexOf('x');
    if (exp_pos == -1) {
        m_degree = 0;
    }
    else {
        if (input.length() <= exp_pos + 1) {
            m_degree = 1;
        }
        else {
            if (input.at(exp_pos + 1) != '^') {
                setError(QObject::QObject::tr("Expected ^, got %1").arg(input.at(exp_pos + 1)));
                return;
            }

            QString exponent = input.right(input.length() - exp_pos - 2);
            m_degree = exponent.toInt(&conversionOk);

            if (!conversionOk) {
                setError(QObject::tr("not valid exponent: %1").arg(exponent));
                return;
            }
        }
    }

    // Parse the coefficient
    parseCoefficient((exp_pos == -1) ? input : input.left(exp_pos));
}

void
Monomial::parseCoefficient(QString coefficient)
{
    if (coefficient.isEmpty()) {
        parseCoefficient("1");
        return;
    }

    // We have to handle a fair number of differente cases in here:
    int indexOfFrac = coefficient.indexOf('/');

    // This case means that we have a fraction, so parse both
    // numerator and denominator.
    if (indexOfFrac != -1) {
        mpq_t t_real, t_imag, zero;
        mpq_init(t_real);
        mpq_init(t_imag);
        mpq_init(zero);

        parseNumber(coefficient.left(indexOfFrac), t_real, t_imag);
        mpq_set(realRationalCoefficient, t_real);
        mpq_set(imagRationalCoefficient, t_imag);

        parseNumber(coefficient.right(coefficient.length() - indexOfFrac - 1), t_real, t_imag);

        // Check that we have not a zero denominator
        if (mpq_equal(t_real, zero) && mpq_equal(t_imag, zero)) {
            setError(QObject::tr("zero denominator in coefficient: %1").arg(coefficient));
            return;
        }

        // Perform complex division
        {
            mpq_t t1, t2;
            mpq_init(t1);
            mpq_init(t2);

            mpq_neg(t_imag, t_imag);
            mpq_mul(t1, t_real, t_real);
            mpq_mul(t2, t_imag, t_imag);

            mpq_add(t1, t1, t2);
            mpq_inv(t1, t1);

            mpq_mul(realRationalCoefficient, realRationalCoefficient, t1);
            mpq_mul(imagRationalCoefficient, imagRationalCoefficient, t1);

            mpq_mul(t1, realRationalCoefficient, t_real);
            mpq_mul(t2, imagRationalCoefficient, t_imag);
            mpq_sub(t1, t1, t2);

            mpq_mul(t2, realRationalCoefficient, t_imag);
            mpq_mul(imagRationalCoefficient, imagRationalCoefficient, t_real);
            mpq_add(imagRationalCoefficient, imagRationalCoefficient, t2);

            mpq_set(realRationalCoefficient, t1);

            mpq_clear(t1);
            mpq_clear(t2);
        }

        mpq_clear(t_real);
        mpq_clear(t_imag);
        mpq_clear(zero);
    }
    else {
        parseNumber(coefficient, realRationalCoefficient, imagRationalCoefficient);
    }
}

void
Monomial::parseNumber(QString number, mpq_t real_output, mpq_t imag_output)
{
    // Sanity check. Check that only allowed symbols are in the number
    foreach (QChar c, number) {
        if (!m_validChars.contains(c)) {
            setError(QObject::tr("Invalid char found: %1").arg(c));
            return;
        }
    }

    // Base case of empty string
    if (number.isEmpty()) {
        mpq_set_str(real_output, "1", 10);
        mpq_set_str(imag_output, "0", 10);
    }

    // Check for the complex case
    if (number.startsWith('(')) {
        // This means that we two numbers to parse, in here.
        // The input may be a string of the type (C + Di) or
        // (C - Di) or ...

        // Start to reduce by +
        int plus_position = number.indexOf('+');
        int minus_position = number.indexOf('-');
        int divider_position = (plus_position != -1) ? plus_position : minus_position;
        if (divider_position != -1) {
            QString part1 = number.left(divider_position) + ')';
            QString part2 = '(' + number.right(number.length() - divider_position - 1);

            mpq_t t_real, t_imag;
            mpq_init(t_real);
            mpq_init(t_imag);

            parseNumber(part1, real_output, imag_output);
            parseNumber(part2, t_real, t_imag);

            if (plus_position != -1) {
                mpq_add(real_output, real_output, t_real);
                mpq_add(imag_output, imag_output, t_imag);
            }
            else {
                mpq_sub(real_output, real_output, t_real);
                mpq_sub(imag_output, imag_output, t_imag);

            }

            mpq_clear(t_real);
            mpq_clear(t_imag);
        }
        else {
            parseNumber(number.mid(1, number.length() - 2), real_output, imag_output);
        }

        return;
    }

    // Check if we are in the pure imaginary case
    if (number.endsWith('i') || number.endsWith('I')) {
        number = number.left(number.length() - 1);
        parseNumber(number, real_output, imag_output);

        // Multiply the number by i
        mpq_t t;
        mpq_init(t);
        mpq_set(t, real_output);
        mpq_neg(real_output, imag_output);
        mpq_set(imag_output, t);
        mpq_clear(t);

        return;
    }

    bool conversionOk;
    int exponent = 0;
    int indexOfE = number.indexOf('e');
    if (indexOfE != -1) {
        // Get the exponent and cut the string
        exponent = number.right(number.length() - indexOfE - 1).toInt(&conversionOk);

        if (!conversionOk) {
            setError(QObject::tr("unable to parse the exponent of the number: %1").arg(number));
            return;
        }

        number = number.left(indexOfE);
    }

    // Now we have a number, that may be composed as a decimal number
    int indexOfDot = number.indexOf('.');
    if (indexOfDot != -1) {
        mpq_t real, imag, t1, t2;
        mpq_init(real);
        mpq_init(imag);
        mpq_init(t1);
        mpq_init(t2);

        QString first = number.left(indexOfDot);
        QString decimal = number.right(number.length() - indexOfDot - 1);

        // Get the leading part first.
        parseNumber(first, real_output, imag_output);

        // ...and then the decimal one. Remeber to divide
        // it for the appropriate number before adding
        // to the leading one.
        parseNumber(decimal, real, imag);
        int exp = decimal.length();

        mpq_set_str(t1, "1", 10);
        mpq_set_str(t2, "10", 10);
        for (int i = 0; i < exp; i++) {
            mpq_mul(t1,t1,t2);
        }

        mpq_div(real, real, t1);
        mpq_div(imag, imag, t1);

        mpq_add(real_output, real_output, real);
        mpq_add(imag_output, imag_output, imag);

        mpq_clear(real);
        mpq_clear(imag);
        mpq_clear(t1);
        mpq_clear(t2);
    }
    else {
        QByteArray data = number.toLocal8Bit();
        mpq_set_str(imag_output, "0", 10);
        mpq_set_str(real_output, data.data(), 10);
    }

    if (exponent != 0) {
        // We have to multiply the number for its exponent.
        mpq_t t1, t2;
        mpq_init(t1);
        mpq_init(t2);

        mpq_set_str(t1, "1", 10);
        mpq_set_str(t2, "10", 10);

        for(int i = 0; i < exponent; i++) {
            mpq_mul(t1, t1, t2);
        }

        mpq_mul(real_output, real_output, t1);
        mpq_mul(imag_output, imag_output, t1);

        mpq_clear(t1);
        mpq_clear(t2);
    }

    mpq_canonicalize(real_output);
    mpq_canonicalize(imag_output);
}

void
Monomial::setError(QString message)
{
    m_errorMessage = message;
    m_valid = FALSE;
}


Monomial::~Monomial()
{
    mpq_clear(realRationalCoefficient);
    mpq_clear(imagRationalCoefficient);
}

int
Monomial::degree() const
{
    return m_degree;
}

QString
Monomial::errorMessage()
{
    return m_errorMessage;
}

bool
Monomial::isValid()
{
    return m_valid;
}

void
Monomial::changeSign()
{
    mpq_neg(realRationalCoefficient, realRationalCoefficient);
    mpq_neg(imagRationalCoefficient, imagRationalCoefficient);
}

void
Monomial::addToMonomialPoly(mps_context *ctx, mps_monomial_poly *poly)
{
    mps_monomial_poly_set_coefficient_q(ctx, poly, degree(),
        realRationalCoefficient, imagRationalCoefficient);
}

Monomial&
Monomial::operator=(const Monomial& rhs)
{
    mpq_set(realRationalCoefficient, rhs.realRationalCoefficient);
    mpq_set(imagRationalCoefficient, rhs.imagRationalCoefficient);

    m_degree = rhs.degree();

    return *this;
}



} // namespace xmpsolve
