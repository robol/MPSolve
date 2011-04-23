/***********************************************************
**            Link library for MT, MPC and GMP            **
**                      Version 1.0                       **
**                                                        **
**             Written by Giuseppe Fiorentino             **
**                 (fiorent@dm.unipi.it)                  **
***********************************************************/

#ifndef __LINK_H__
#define __LINK_H__

/* needed header files */
#include <mps/mt.h>
#include <mps/mpc.h>

/***********************************************************
**              link functions                            **
***********************************************************/

void mpf_set_rdpe(mpf_t f, rdpe_t e);
void mpf_get_rdpe(rdpe_t e, mpf_t f);

void mpc_set_cplx(mpc_t mc, cplx_t c);
void mpc_get_cplx(cplx_t c, mpc_t mc);

void mpc_set_cdpe(mpc_t mc, cdpe_t c);
void mpc_get_cdpe(cdpe_t c, mpc_t mc);

#endif

/***********************************************************
**                                                        **
***********************************************************/
