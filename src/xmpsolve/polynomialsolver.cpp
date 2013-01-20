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
    PolynomialParser parser(m_mpsContext);

    // Parse the input string that the user has given.
    mps_monomial_poly * poly = parser.parse(inputString);

    if (poly != NULL) {
        return solvePoly(poly);
    }
    else {
        qDebug() << tr("Cannot parse input polynomial string: %1").arg(parser.errorMessage());
        return -1;
    }
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
