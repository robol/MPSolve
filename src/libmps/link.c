/***********************************************************
**            Link library for MT, MPC and GMP            **
**                      Version 1.1                       **
**                                                        **
**             Written by Giuseppe Fiorentino             **
**                 (fiorent@dm.unipi.it)                  **
***********************************************************/

#include <mps/link.h>

/***********************************************************
**              link functions                            **
***********************************************************/

void
mpf_set_rdpe (mpf_t f, rdpe_t e)
{
  mpf_set_d (f, rdpe_Mnt (e));
  if (rdpe_Esp (e) >= 0)
    mpf_mul_2exp (f, f, rdpe_Esp (e));
  else
    mpf_div_2exp (f, f, -rdpe_Esp (e));
}

void
mpf_get_rdpe (rdpe_t e, mpf_t f)
{
  mp_exp_t esp;

  esp = f->_mp_exp;
  f->_mp_exp = 0;
  rdpe_set_2dl (e, mpf_get_d (f), esp * mp_bits_per_limb);
  f->_mp_exp = esp;
}

void
mpc_set_cplx (mpc_t mc, cplx_t c)
{
  mpf_set_d (mpc_Re (mc), cplx_Re (c));
  mpf_set_d (mpc_Im (mc), cplx_Im (c));
}

void
mpc_get_cplx (cplx_t c, mpc_t mc)
{
  cplx_Re (c) = mpf_get_d (mpc_Re (mc));
  cplx_Im (c) = mpf_get_d (mpc_Im (mc));
}

void
mpc_set_cdpe (mpc_t mc, cdpe_t c)
{
  mpf_set_rdpe (mpc_Re (mc), cdpe_Re (c));
  mpf_set_rdpe (mpc_Im (mc), cdpe_Im (c));
}

void
mpc_get_cdpe (cdpe_t c, mpc_t mc)
{
  mpf_get_rdpe (cdpe_Re (c), mpc_Re (mc));
  mpf_get_rdpe (cdpe_Im (c), mpc_Im (mc));
}

/***********************************************************
**                                                        **
***********************************************************/
