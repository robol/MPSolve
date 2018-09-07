#include "polynomialsolver.h"
#include "mpsolveworker.h"
#include "root.h"
#include <QDebug>
#include <QRegExp>
#include <QStringList>

#include <QDebug>

static double LOG2_10 = log(10) / log(2);

using namespace xmpsolve;

PolynomialSolver::PolynomialSolver(QObject *parent) :
    QObject(parent)
{
    m_mpsContext = NULL;
    m_currentPoly = NULL;
    m_worker.connect(&m_worker, SIGNAL(finished()),
                     this, SLOT(workerExited()));
    m_errorMessage = "";
}

PolynomialSolver::~PolynomialSolver()
{
    if (m_currentPoly)
      mps_polynomial_free (m_mpsContext, m_currentPoly);
    if (m_mpsContext)
      mps_context_free(m_mpsContext);
}

int
PolynomialSolver::solvePolFile(QString selectedFile, mps_algorithm selected_algorithm, int required_digits,
                               mps_output_goal goal)
{
    QByteArray stringData = selectedFile.toLatin1().data();

    m_mpsContext = mps_context_new();
    m_worker.setMpsContext(m_mpsContext);

    // Parse the stream specified by the user
    mps_polynomial * poly = mps_parse_file (m_mpsContext, stringData.data());

    if (mps_context_has_errors (m_mpsContext) || !poly) {
        m_errorMessage = tr("Error while solving the given pol file: %1").
                arg(mps_context_error_msg(m_mpsContext));
        mps_context_free (m_mpsContext);
        m_mpsContext = NULL;
        return -1;
    }
    else
        mps_context_set_input_poly (m_mpsContext, poly);

    // Select the options selected by the user
    mps_context_select_algorithm(m_mpsContext, selected_algorithm);
    mps_context_set_output_prec(m_mpsContext, required_digits * LOG2_10);
    mps_context_set_output_goal(m_mpsContext, goal);

    m_worker.start();
    return mps_context_get_degree (m_mpsContext);
}

int
PolynomialSolver::solvePolFileFromContent(QString content, mps_algorithm selected_algorithm, int required_digits,
                                          mps_output_goal goal)
{
    QByteArray stringData = content.toLatin1().data();

    m_mpsContext = mps_context_new();
    m_worker.setMpsContext(m_mpsContext);

    // Parse the stream specified by the user
    mps_polynomial * poly = mps_parse_string (m_mpsContext, stringData.data());

    if (mps_context_has_errors (m_mpsContext) || !poly) {
        m_errorMessage = tr("Error while solving the given pol file: %1").
                arg(mps_context_error_msg(m_mpsContext));
        mps_context_free (m_mpsContext);
        m_mpsContext = NULL;
        return -1;
    }
    else
        mps_context_set_input_poly (m_mpsContext, poly);

    // Select the options selected by the user
    mps_context_select_algorithm(m_mpsContext, selected_algorithm);
    mps_context_set_output_prec(m_mpsContext, required_digits * LOG2_10);
    mps_context_set_output_goal(m_mpsContext, goal);

    m_worker.start();
    return mps_context_get_degree (m_mpsContext);
}

int
PolynomialSolver::solvePoly(mps_polynomial * poly, PolynomialBasis basis,
                            mps_algorithm selected_algorithm,
                            int required_digits, mps_output_goal goal)
{
    m_currentPoly = poly;
    
    mps_context_set_input_poly(m_mpsContext, MPS_POLYNOMIAL (m_currentPoly));
    mps_context_select_algorithm(m_mpsContext, selected_algorithm);
    mps_context_set_output_prec(m_mpsContext, required_digits * LOG2_10);
    mps_context_set_output_goal(m_mpsContext, goal);

    // One might want to uncomment this for debugging purposes. 
    // mps_context_add_debug_domain (m_mpsContext, MPS_DEBUG_TRACE);
    
    m_worker.setMpsContext (m_mpsContext);

    m_worker.start();
    return mps_context_get_degree (m_mpsContext);
}

int
PolynomialSolver::solvePoly(QString inputString, PolynomialBasis basis,
                            mps_algorithm selected_algorithm,
                            int required_digits, mps_output_goal goal)
{
  if (! m_mpsContext)
    m_mpsContext = mps_context_new();

  QByteArray inputStringData = inputString.toLatin1();
  mps_polynomial * poly = (mps_polynomial *) mps::Polynomial::fromString (m_mpsContext, inputStringData.data());
  
  // If the user has chosen the Chebyshev basis, keep the coefficients but
  // change the basis of the parsed polynomial. 
  if (basis == PolynomialBasis::CHEBYSHEV) {
	mpq_t real_coeff, imag_coeff;
	mps_monomial_poly *mpoly = MPS_MONOMIAL_POLY (poly);
	  
	mps_chebyshev_poly* cpoly = mps_chebyshev_poly_new (m_mpsContext, 
	  poly->degree, poly->structure);
	  
	// Copy the coefficients -- every polynomial parsed from a command line has
	// rational coefficients. 
	mpq_init (real_coeff);
	mpq_init (imag_coeff);
	
	for (int i = 0; i <= poly->degree; i++) {
	  mps_monomial_poly_get_coefficient_q  (m_mpsContext, mpoly, i, real_coeff, imag_coeff);
	  mps_chebyshev_poly_set_coefficient_q (m_mpsContext, cpoly, i, real_coeff, imag_coeff);
	}
	
	mpq_clear (real_coeff);
	mpq_clear (imag_coeff);
	
	poly = MPS_POLYNOMIAL (cpoly);
	mps_polynomial_free (m_mpsContext, MPS_POLYNOMIAL (mpoly));
  }

  if (poly && (! mps_context_has_errors(m_mpsContext)) && poly->degree != 0) {
    return solvePoly(poly, basis, selected_algorithm, required_digits, goal);
  }
  else {
    mps_context_free (m_mpsContext);
    m_mpsContext = NULL;
    
    m_errorMessage = tr("There is a syntax error in the specified the polynomial");
    
    if (poly && poly->degree == 0) {
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
    mpc_t * results = NULL;
    rdpe_t * radii  = NULL;

    mps_context_get_roots_m(m_mpsContext, &results, &radii);

    for (int i = 0; i < n; i++)
    {
        mps_root_status status = mps_context_get_root_status (m_mpsContext, i);
        Root* r = new Root(results[i], radii[i], status);
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

    if (results != NULL)
      free (results);
    if (radii != NULL)
      free (radii);

    // Update the internal model with the approximations
    m_rootsModel.setRoots(roots);

    /* Set the roots somewhere and then called
     * the solved() signal */
    solved();
}

void
PolynomialSolver::abortComputations()
{
    m_worker.abortComputation();
}

unsigned long int
PolynomialSolver::CPUTime()
{
    return m_worker.CPUTime();
}

RootsModel *PolynomialSolver::rootsModel()
{
    return &m_rootsModel;
}
