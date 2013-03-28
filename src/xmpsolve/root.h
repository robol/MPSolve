#ifndef xmpsolve_ROOT_H
#define xmpsolve_ROOT_H

#include <mps/mps.h>

namespace xmpsolve {

class Root
{

public:
    explicit Root();
    explicit Root(mpc_t value, rdpe_t radius, mps_root_status status = MPS_ROOT_STATUS_CLUSTERED);
    double get_real_part();
    double get_imag_part();
    double get_radius();

    mpc_t value;
    rdpe_t radius;
    mps_root_status status;
};

} // namespace xmpsolve

#endif // xmpsolve_ROOT_H
