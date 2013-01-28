#include "mpsolveworker.h"
#include <mps/mps.h>
#include <QDebug>

using namespace xmpsolve;

MPSolveWorker::MPSolveWorker(mps_context * s, QObject *parent) :
    QThread(parent)
{
    /* Set the status in the worker */
    m_context = s;
}

void
MPSolveWorker::run()
{
    /* Actually solve the polynomial that should have been
     * set in here... */
    clock_t *timer = mps_start_timer();
    mps_mpsolve(m_context);
    m_time = mps_stop_timer(timer);
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
