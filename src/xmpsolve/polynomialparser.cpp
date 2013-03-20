#include "polynomialparser.h"
#include "polynomial.h"
#include <QDebug>

namespace xmpsolve {

PolynomialParser::PolynomialParser(QObject *parent) :
    QObject(parent)
{
    m_degree = 0;
}

void
PolynomialParser::reset()
{
    m_degree = 0;
    qDeleteAll(m_monomials);
    m_monomials.clear();
}

Polynomial
PolynomialParser::parse(QString input)
{
    bool signChangeNeeded = false;

    if (input.isEmpty()) {
        m_errorMessage = tr("Empty string encountered while parsing.");
        return Polynomial();
    }

    input = input.trimmed();

    if (input.at(0) == '-') {
        signChangeNeeded = true;
        input = input.right(input.length() - 1).trimmed();
    }
    else if (input.at(0) == '+') {
        input = input.right(input.length() - 1).trimmed();
    }

    // We need to tokenize the input string dividing each monomial.
    // We do that scanning recursively for them, tokenizing by + and -
    int position = 0;
    int bracket_depth = 0;
    while (position < input.length() && (
           (input.at(position) != '+' &&
           input.at(position) != '-') ||
           bracket_depth != 0))
    {
        if (input.at(position) == '(') {
            bracket_depth++;
        }
        if (input.at(position) == ')') {
            bracket_depth--;
        }

        position++;
    }

    if (bracket_depth != 0) {
        m_errorMessage = tr("Mismatched parenthesis");
        return Polynomial();
    }

    // So now we have that input[position] \in { +, - } or
    // position == input.length(). In both cases we parse the
    // monomial described by input[:position].
    Monomial *monomial = new Monomial(input.left(position));

    if (!monomial->isValid()) {
        m_errorMessage = monomial->errorMessage();
        return Polynomial();
    }

    if (signChangeNeeded) {
        monomial->changeSign();
    }

    if (monomial->degree() > m_degree)
        m_degree = monomial->degree();

    m_monomials.append(monomial);

    if (position < input.length()) {
        return parse(input.right(input.length() - position));
    }
    else {
        Polynomial sumPoly;

        foreach(Monomial *m, m_monomials) {
            sumPoly += *m;
        }

        return sumPoly;
    }
}

QString
PolynomialParser::errorMessage()
{
    return m_errorMessage;
}

} // namespace xmpsolve
