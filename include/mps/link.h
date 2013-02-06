/***********************************************************
**            Link library for MT, MPC and GMP            **
**                      Version 1.0                       **
**                                                        **
**             Written by Giuseppe Fiorentino             **
**                 (fiorent@dm.unipi.it)                  **
***********************************************************/

/**
 * @file
 * @brief Functions used to pass the data back and forth between
 * multiprecision and floating point.
 */

#ifndef __LINK_H__
#define __LINK_H__

/* needed header files */
#include <mps/mt.h>
#include <mps/mpc.h>

#ifdef __cplusplus
extern "C" {
#endif

/***********************************************************
**              link functions                            **
***********************************************************/

/**
 * @brief Set the Multiprecision value <code>f</code> with the value
 * stored in <code>e</code>.
 *
 * @param f The multiprecision floating point number to set.
 * @param e The RDPE value to set in <code>f</code>.
 */
void mpf_set_rdpe (mpf_t f, rdpe_t e);

/**
 * @brief Get the RDPE version of the Multiprecision value <code>f</code>.
 *
 * @param e The RDPE where the value of <code>f</code> will be stored.
 * @param f The multiprecision floating point number to extract the value from.
 */
void mpf_get_rdpe (rdpe_t e, mpf_t f);

/**
 * @brief Set the Multiprecision value <code>f</code> with the value
 * stored in <code>e</code>.
 *
 * @param mc The multiprecision complex number to set.
 * @param c The <code>cplx_t</code> value to set in <code>f</code>.
 */
void mpc_set_cplx (mpc_t mc, cplx_t c);

/**
 * @brief Get the <code>cplx_t</code> version of the Multiprecision value <code>f</code>.
 *
 * @param c The <code>cplx_t</code> where the value of <code>f</code> will be stored.
 * @param mc The multiprecision complex number to extract the value from.
 */
void mpc_get_cplx (cplx_t c, mpc_t mc);

/**
 * @brief Set the Multiprecision value <code>f</code> with the value
 * stored in <code>e</code>.
 *
 * @param mc The multiprecision complex number to set.
 * @param c The CDPE value to set in <code>f</code>.
 */
void mpc_set_cdpe (mpc_t mc, cdpe_t c);

/**
 * @brief Get the CDPE version of the Multiprecision value <code>f</code>.
 *
 * @param c The CDPE where the value of <code>f</code> will be stored.
 * @param mc The multiprecision complex number to extract the value from.
 */
void mpc_get_cdpe (cdpe_t c, mpc_t mc);

#ifdef __cplusplus
}
#endif

#endif

/***********************************************************
**                                                        **
***********************************************************/
