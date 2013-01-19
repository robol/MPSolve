#ifndef MPSOLVEWORKER_H
#define MPSOLVEWORKER_H

#include <QThread>
#include <mps/mps.h>

namespace xmpsolve {

class MPSolveWorker : public QThread
{
    Q_OBJECT
public:
    explicit MPSolveWorker(mps_context * s, QObject *parent = 0);
    void run();

private:
    mps_context * m_context;
    
signals:
    
public slots:
    
};

} // Namespace xmpsolve

#endif // MPSOLVEWORKER_H
