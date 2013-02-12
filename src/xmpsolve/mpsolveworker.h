#ifndef MPSOLVEWORKER_H
#define MPSOLVEWORKER_H

#include <QThread>
#include <mps/mps.h>

namespace xmpsolve {

class MPSolveWorker : public QThread
{
    Q_OBJECT
public:
    explicit MPSolveWorker(mps_context * s = NULL, QObject *parent = 0);

    /**
     * @brief setMpsContext can be used to set the current mps_context
     * for the computation.
     * @param ctx is the pointer to the current mps_context.
     */
    void setMpsContext(mps_context * ctx);

    /**
     * @brief run Actually start the computation calling mps_mpsolve();
     */
    void run();

    /**
     * @brief CPUTime gets the number of ms of CPU time used by the
     * last call to run().
     * @return The number of ms of CPU time spent on the last run.
     */
    unsigned long int CPUTime();

private:
    mps_context * m_context;
    unsigned long int m_time;
    
signals:
    
public slots:
    
};

} // Namespace xmpsolve

#endif // MPSOLVEWORKER_H
