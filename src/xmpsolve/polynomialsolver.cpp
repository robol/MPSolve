#include "polynomialsolver.h"
#include "mpsolveworker.h"
#include "root.h"
#include <QDebug>
#include <QRegExp>
#include <QStringList>

using namespace xmpsolve;

PolynomialSolver::PolynomialSolver(QObject *parent) :
    QObject(parent)
{
    m_mpsContext = mps_context_new();
}

PolynomialSolver::~PolynomialSolver()
{
    mps_context_free(m_mpsContext);
}

int
PolynomialSolver::solvePoly(mps_monomial_poly *poly)
{
    mps_context_set_input_poly(m_mpsContext, poly);
    mps_context_select_algorithm(m_mpsContext, MPS_ALGORITHM_SECULAR_GA);
    mps_context_set_output_goal(m_mpsContext, MPS_OUTPUT_GOAL_APPROXIMATE);

    // Create a new thread that solve the polynomial.
    m_worker = new MPSolveWorker(m_mpsContext);
    m_worker->connect(m_worker, SIGNAL(finished()),
                      this, SLOT(workerExited()));

    m_worker->start();
    return mps_context_get_degree (m_mpsContext);
}

int
PolynomialSolver::solvePoly(QString inputString)
{
    int degree = -1;

    if (!inputString.isEmpty())
        degree = 0;
    else
        return degree;

    /* Find out the degree of the polynomial */
    QRegExp exponents("x\\^(\\d+)");

    bool conversion_ok;
    int pos = 0;
    while ((pos = exponents.indexIn(inputString, pos)) != -1)
    {
        int new_degree = exponents.cap(1).toInt(&conversion_ok);
        if (!conversion_ok)
        {
           degree = -1;
           break;
        }
        if (new_degree > degree)
            degree = new_degree;
        pos += exponents.matchedLength();
    }

    // qDebug() << "Degree = " << degree;

    /* Now allocate the polynomial and find the coefficients */
    mps_monomial_poly * p = mps_monomial_poly_new(m_mpsContext, degree);

    int sign = 1;

    /* Scan the input for coefficients */
    int aux_pos;
    for (pos = 0; pos < inputString.length(); pos++)
    {

        if (inputString[pos].isSpace())
        {
            continue;
        }

        if (inputString[pos] == '-')
        {
            sign = -1;
            continue;
        }

        /* Start of a coefficient */
        if (inputString[pos].isNumber() || inputString[pos] == 'x' || inputString[pos] == 'X')
        {
            QString coeff;
            bool imag = false;
            aux_pos = pos;

            if (inputString[pos] == 'x' || inputString[pos] == 'X')
            {
                coeff = "1";
            }
            else
            {

                while ((inputString[aux_pos].isNumber() || inputString[aux_pos] == '.' || inputString[aux_pos] == 'e'
                        || inputString[aux_pos] == 'E' || inputString[aux_pos] == 'i') && aux_pos < inputString.length())
                {
                    aux_pos++;
                    if (inputString[aux_pos] == 'i')
                        imag = true;
                }

                coeff = inputString.mid (pos, aux_pos - pos);
            }

            /* Find out the degree */
            QString d;
            int c_pos;
            aux_pos++;
            while (!(inputString[aux_pos] == '^' || aux_pos >= inputString.length() ||
                   inputString[aux_pos].isSpace() || inputString[aux_pos] == '-'
                   || inputString[aux_pos] == '+'))
                aux_pos++;
            aux_pos++;

            c_pos = aux_pos;

            if (inputString[aux_pos - 1].isSpace() || aux_pos - 1 == inputString.length()
                || inputString[aux_pos - 1] == '-' || inputString[aux_pos - 1]  == '+')
            {
                if (aux_pos >= 2 && inputString[aux_pos - 2] == 'x')
                    d = "1";
                else
                    d = "0";
            }
            else
            {

                while (inputString[c_pos].isNumber() && c_pos <= inputString.length())
                {
                    c_pos++;
                }
                c_pos++;

                d = inputString.mid (aux_pos, c_pos - aux_pos);
            }

            int new_degree = d.toInt(&conversion_ok);
            mps_monomial_poly_set_coefficient_d(m_mpsContext, p, new_degree,
                                                sign * coeff.toFloat(), 0);

            pos = c_pos - 2;
            sign = +1;
        }
    }

    return solvePoly(p);
}

void
PolynomialSolver::workerExited()
{
    int n = mps_context_get_degree (m_mpsContext);

    /* Save the roots in a list */
    QList<Root*> roots;

    /* Prepare some space to save the roots */
    mpc_t * results = new mpc_t[n];
    rdpe_t * radii = new rdpe_t[n];

    mpc_vinit2 (results, n, 53);

    mps_context_get_roots_m(m_mpsContext, results, radii);

    for (int i = 0; i < n; i++)
    {
        Root* r = new Root();

        mpc_init2 (r->value, mpc_get_prec (results[i]));
        mpc_set (r->value, results[i]);
        rdpe_set (r->radius, radii[i]);

        roots.append(r);
    }

    mpc_vclear (results, n);
    delete [] results;
    delete [] radii;

    /* Set the roots somewhere and then called
     * the solved() signal */
    solved(roots);
}
