#include "root.h"

namespace xmpsolve {

Root::Root()
{
    mpcf_init2  (this->value, 0);
    mpcf_set_ui (this->value, 0U, 0U);

    rdpe_set (this->radius, rdpe_zero);
    this->status = MPS_ROOT_STATUS_CLUSTERED;
}

Root::Root(mpcf_t value, rdpe_t radius, mps_root_status status)
{
    mpcf_init2 (this->value, mpcf_get_prec (value));
    mpcf_set  (this->value, value);

    rdpe_set (this->radius, radius);

    this->status = status;
}

double
Root::get_radius()
{
    return rdpe_get_d(radius);
}

double
Root::get_real_part()
{
    return mpf_get_d (mpcf_Re (value));
}

double
Root::get_imag_part()
{
    return mpf_get_d (mpcf_Im (value));
}

} // namespace xmpsolve
