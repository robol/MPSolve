/*
 * This file is part of MPSolve 3.2.1
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */


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
mpcf_set_cplx (mpcf_t mc, cplx_t c)
{
  mpf_set_d (mpcf_Re (mc), cplx_Re (c));
  mpf_set_d (mpcf_Im (mc), cplx_Im (c));
}

void
mpcf_get_cplx (cplx_t c, mpcf_t mc)
{
#ifdef MPS_USE_BUILTIN_COMPLEX
  cplx_Re (c) = mpf_get_d (mpcf_Re (mc));
  cplx_Im (c) = mpf_get_d (mpcf_Im (mc));
#else
  cplx_set_d (c, mpf_get_d (mpcf_Re (mc)),
              mpf_get_d (mpcf_Im (mc)));
#endif
}

void
mpcf_set_cdpe (mpcf_t mc, cdpe_t c)
{
  mpf_set_rdpe (mpcf_Re (mc), cdpe_Re (c));
  mpf_set_rdpe (mpcf_Im (mc), cdpe_Im (c));
}

void
mpcf_get_cdpe (cdpe_t c, mpcf_t mc)
{
  mpf_get_rdpe (cdpe_Re (c), mpcf_Re (mc));
  mpf_get_rdpe (cdpe_Im (c), mpcf_Im (mc));
}

/***********************************************************
**                                                        **
***********************************************************/
