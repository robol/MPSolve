/* 
 * File:   secular.h
 * Author: leonardo
 *
 * Created on 29 aprile 2011, 16.26
 */

#ifndef SECULAR_H
#define	 SECULAR_H

#ifdef	__cplusplus
extern "C"
{
#endif

#include <mps/mt.h>
#include <mps/interface.h>
#include <mps/mpc.h>
#include <float.h>

/* MACROS */
#define mps_secular_equation_from_status(s) (mps_secular_equation*) (s)->secular_equation

/* Routines in secular-newton.c */
  void
    mps_secular_fnewton (mps_status * st, cplx_t x, double *rad, cplx_t corr,
			 mps_boolean * again);
  void
    mps_secular_dnewton (mps_status * st, cdpe_t x, rdpe_t rad, cdpe_t corr,
			 mps_boolean * again);
  void
    mps_secular_mnewton (mps_status * st, mpc_t x, rdpe_t rad, mpc_t corr,
			 mps_boolean * again);

/* Routines in secular.c */
  void mps_secular_deflate (mps_status * s, mps_secular_equation * sec);

  void mps_secular_check_data (mps_status * s, char *which_case);

  void mps_secular_switch_phase (mps_status * s, mps_phase phase);

  void mps_secular_raise_coefficient_precision (mps_status * s, int wp);

  void mps_secular_raise_precision (mps_status * s, int wp);

/* Routines in secular-starting.c */
  void
    mps_secular_fstart (mps_status * s, int n, int i_clust, double clust_rad,
			double g, rdpe_t eps);
  void
    mps_secular_dstart (mps_status * s, int n, int i_clust, rdpe_t clust_rad,
			rdpe_t g, rdpe_t eps);
  void
    mps_secular_mstart (mps_status * s, int n, int i_clust, rdpe_t clust_rad,
			rdpe_t g, rdpe_t eps);

/* Routines in secular-ga.c */
  int mps_secular_ga_fiterate (mps_status * s, int maxit);

  int mps_secular_ga_diterate (mps_status * s, int maxit);

  int mps_secular_ga_miterate (mps_status * s, int maxit);

  void mps_secular_ga_regenerate_coefficients (mps_status * s);

  mps_boolean mps_secular_ga_check_stop (mps_status * s);

  void mps_secular_ga_mpsolve (mps_status * s);

/* Interface functions in mps_secular.c */
  mps_secular_equation *mps_secular_equation_new (mps_status * s,
						  cplx_t * afpc,
						  cplx_t * bfpc,
						  unsigned long int n);

  mps_secular_equation *mps_secular_equation_new_raw (mps_status * s,
						      unsigned long int n);

  void mps_secular_equation_free (mps_secular_equation * s);

  void mps_secular_set_radii (mps_status * s);

#ifdef	__cplusplus
}
#endif

#endif				/* SECULAR_H */
