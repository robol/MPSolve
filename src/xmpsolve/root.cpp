#include "root.h"

namespace xmpsolve {

Root::Root(QObject *parent) :
    QObject(parent)
{
}

double
Root::get_radius()
{
    return rdpe_get_d(radius);
}

double
Root::get_real_part()
{
    return mpf_get_d (mpc_Re (value));
}

double
Root::get_imag_part()
{
    return mpf_get_d (mpc_Im (value));
}

} // namespace xmpsolve
