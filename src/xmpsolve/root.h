#ifndef xmpsolve_ROOT_H
#define xmpsolve_ROOT_H

#include <QObject>
#include <mps/mps.h>

namespace xmpsolve {

class Root : public QObject
{
    Q_OBJECT
public:
    explicit Root(QObject *parent = 0);
    double get_real_part();
    double get_imag_part();
    double get_radius();

    mpc_t value;
    rdpe_t radius;
    
signals:
    
public slots:
    
};

} // namespace xmpsolve

#endif // xmpsolve_ROOT_H
