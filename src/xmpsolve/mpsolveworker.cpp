#include "mpsolveworker.h"
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
    mps_mpsolve(m_context);
}
