#include "polynomialsolver.h"
#include "mpsolveworker.h"
#include "root.h"
#include "polynomialparser.h"
#include <QDebug>
#include <QRegExp>
#include <QStringList>

using namespace xmpsolve;

PolynomialSolver::PolynomialSolver(QObject *parent) :
    QObject(parent)
{
    m_mpsContext = NULL;
    m_worker.connect(&m_worker, SIGNAL(finished()),
                     this, SLOT(workerExited()));
}

PolynomialSolver::~PolynomialSolver()
{
    if (m_mpsContext)
        mps_context_free(m_mpsContext);
}

int
PolynomialSolver::solvePoly(Polynomial poly, mps_algorithm selected_algorithm,
                            int required_digits)
{
    m_currentPoly = poly;

    // Regenerate a new mps_context, since as of now the same cannot be used
    // more than one time.
    if (m_mpsContext)
        mps_context_free(m_mpsContext);
    m_mpsContext = mps_context_new();

    m_worker.setMpsContext(m_mpsContext);

    mps_monomial_poly * monomialPoly = mps_monomial_poly_new(m_mpsContext, poly.degree());
    for (int i = 0; i <= poly.degree(); i++) {
        poly.monomial(i).addToMonomialPoly(m_mpsContext, monomialPoly);
    }

    mps_context_set_input_poly(m_mpsContext, MPS_POLYNOMIAL (monomialPoly));
    mps_context_select_algorithm(m_mpsContext, selected_algorithm);
    mps_context_set_output_prec(m_mpsContext, required_digits * LOG2_10);
    mps_context_set_output_goal(m_mpsContext, MPS_OUTPUT_GOAL_APPROXIMATE);

    m_worker.start();
    return mps_context_get_degree (m_mpsContext);
}

int
PolynomialSolver::solvePoly(QString inputString, mps_algorithm selected_algorithm,
                            int required_digits)
{
    PolynomialParser parser;

    // Parse the input string that the user has given.
    Polynomial poly = parser.parse(inputString);

    if (poly.degree() != 0) {
        return solvePoly(poly, selected_algorithm, required_digits);
    }
    else {
       m_errorMessage = parser.errorMessage();

       if (m_errorMessage == QString("")) {
           m_errorMessage = tr("Constant polynomials have no roots");
       }

       return -1;
    }
}

QString
PolynomialSolver::errorMessage()
{
    return m_errorMessage;
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

    mps_context_get_roots_m(m_mpsContext, &results, &radii);

    for (int i = 0; i < n; i++)
    {
        Root* r = new Root();

        mpc_init2 (r->value, mpc_get_prec (results[i]));
        mpc_set (r->value, results[i]);
        rdpe_set (r->radius, radii[i]);

        roots.append(r);
    }

    for (int i = 0; i < mps_context_get_zero_roots (m_mpsContext); i++)
    {
        Root * r = new Root();

        mpc_init2 (r->value, 0);
        mpc_set_ui (r->value, 0U, 0U);
        rdpe_set (r->radius, rdpe_zero);

        roots.append(r);
    }

    mpc_vclear (results, n);
    delete [] results;
    delete [] radii;

    /* Set the roots somewhere and then called
     * the solved() signal */
    solved(roots);
}

unsigned long int
PolynomialSolver::CPUTime()
{
    return m_worker.CPUTime();
}
