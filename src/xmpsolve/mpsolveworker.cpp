#include "mpsolveworker.h"
#include <mps/mps.h>
#include <QDebug>

using namespace xmpsolve;

MPSolveWorker::MPSolveWorker(mps_context * s, QObject *parent) :
    QThread(parent)
{
    /* Set the status in the worker */
    m_context = s;
    m_time = 0;
}

void
MPSolveWorker::run()
{
    /* Actually solve the polynomial that should have been
     * set in here... */
    m_timer = mps_start_timer();
    mps_mpsolve(m_context);
}

void
MPSolveWorker::abortComputation()
{
    mps_context_abort(m_context);
}

unsigned long int
MPSolveWorker::CPUTime()
{
    return m_time;
}

void
MPSolveWorker::setMpsContext(mps_context *ctx)
{
    m_context = ctx;
}
