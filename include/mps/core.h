/***********************************************************
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.2, May 2001                  **
 **                                                        **
 **                      Written by                        **
 **       Dario Andrea Bini and Giuseppe Fiorentino        **
 **       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
 **                                                        **
 ** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
 ***********************************************************/

/**
 * @file
 *
 * This file is the header for the libmps library. Including
 * this file is needed to access all the MPSolve routines by
 * MPSolve internals.
 *
 * If you need an interface for using mps as a library try to
 * include the header file <code>mps/interface.h</code>.
 *
 * @brief Header file for libmps
 */

#ifndef MPS_CORE_H_
#define MPS_CORE_H_

#ifdef __cplusplus
extern "C"
{
#endif


/* local include files */
#include <mps/types.h>
#include <mps/tools.h>
#include <mps/mt.h>
#include <mps/gmptools.h>
#include <mps/mpc.h>
#include <mps/link.h>
#include <mps/debug.h>

/* Interface should be a subset of core, so what is defined
 * there should be included here. */
#include <mps/interface.h>
#include <mps/threading.h>

/* constants */

#define MPS_ALL_CLUSTERS -1


/* FUNCTIONS */

/* Functions in input-buffer.c */
#define MPS_INPUT_BUFFER_HISTORY_DEFAULT_SIZE 2
  mps_input_buffer *mps_input_buffer_new (FILE * stream);
  int mps_input_buffer_readline (mps_input_buffer * buf);
  void mps_input_buffer_free (mps_input_buffer * buf);
  void mps_input_buffer_set_history_size (mps_input_buffer * buf, size_t size);
  mps_boolean mps_input_buffer_eof (mps_input_buffer * buf);
  char * mps_input_buffer_next_token (mps_input_buffer * buf);

/* functions in aberth.c */
  void mps_faberth (mps_status * s, int j, cplx_t abcorr);
  void mps_daberth (mps_status * s, int j, cdpe_t abcorr);
  void mps_maberth (mps_status * s, int j, mpc_t abcorr);
  void mps_faberth_s (mps_status * s, int j, mps_cluster * cluster, cplx_t abcorr);
  void mps_daberth_s (mps_status * s, int j, mps_cluster * cluster, cdpe_t abcorr);
  void mps_maberth_s (mps_status * s, int j, mps_cluster * cluster, mpc_t abcorr);
  void mps_maberth_s_wl (mps_status * s, int j, mps_cluster * cluster, mpc_t abcorr,
                         pthread_mutex_t * aberth_mutex);
  void mps_mnewtis (mps_status * s);

/* functions in cluster.c */
  void mps_cluster_reset (mps_status * s);
  void mps_fcluster (mps_status * s, double * frad, int nf);
  void mps_dcluster (mps_status * s, int nf);
  void mps_mcluster (mps_status * s, int nf);
  /* void mps_cluster_detach (mps_status * s, mps_cluster * cluster); */
  /* void mps_cluster_reassemble (mps_status * s, mps_cluster * cluster); */
  void mps_debug_cluster_structure (mps_status * s);

/* functions in convex.c */
  void mps_fconvex (mps_status * s, int n, double a[]);

/* functions in data.c */
  void mps_mp_set_prec (mps_status * s, long int prec);
  void mps_allocate_data (mps_status * s);
  void mps_prepare_data (mps_status * s, long int prec);
  void mps_restore_data (mps_status * s);
  void mps_free_data (mps_status * s);
  void mps_raise_data_raw (mps_status * s, long int prec);

/* functions in improve.c */
  void mps_improve (mps_status * s);

/* functions in main.c */
  void mps_setup (mps_status * s);
  void mps_check_data (mps_status * s, char *which_case);
  void mps_compute_sep (mps_status * s);
  void mps_standard_mpsolve (mps_status * s);

/* functions in newton.c */
  void mps_fnewton (mps_status * st, int n, cplx_t z, double *radius,
                    cplx_t corr, cplx_t fpc[], double fap[],
                    mps_boolean * cont, mps_boolean skip_radius_computation);
  void mps_dnewton (mps_status * st, int n, cdpe_t z, rdpe_t radius,
                    cdpe_t corr, cdpe_t dpc[], rdpe_t dap[],
                    mps_boolean * cont, mps_boolean skip_radius_computation);
  void mps_mnewton (mps_status * st, int n, mpc_t z, rdpe_t radius,
                    mpc_t corr, mpc_t mfpc[], mpc_t mfppc[], rdpe_t dap[],
                    mps_boolean * spar, mps_boolean * cont, int n_thread, 
		    mps_boolean skip_radius_computation);
  void mps_parhorner (mps_status * st, int n, mpc_t x, mpc_t p[],
                      mps_boolean b[], mpc_t s, int n_thread);
  void mps_aparhorner (mps_status * st, int n, rdpe_t x, rdpe_t p[],
                       mps_boolean b[], rdpe_t s, int n_thread);

  /* Functions in general-radius.c */
  void mps_fradii (mps_status * s, double * fradii);
  void mps_dradii (mps_status * s, rdpe_t * dradii);
  void mps_mradii (mps_status * s, rdpe_t * dradii);

  /* Functions in monomial-radius.c */
  void mps_monomial_fradii (mps_status * s, double * fradii);
  void mps_monomial_dradii (mps_status * s, rdpe_t * dradii);
  void mps_monomial_mradii (mps_status * s, rdpe_t * dradii);

  /* Functions in secular-radius.c */
  void mps_secular_fradii (mps_status * s, double * fradii);
  void mps_secular_dradii (mps_status * s, rdpe_t * dradii);
  void mps_secular_mradii (mps_status * s, rdpe_t * dradii);

  /* Functions in secular-evaluation.c */
  void mps_secular_feval (mps_status * s, mps_secular_equation * sec, cplx_t x, cplx_t value);
  void mps_secular_deval (mps_status * s, mps_secular_equation * sec, cdpe_t x, cdpe_t value);
  void mps_secular_meval (mps_status * s, mps_secular_equation * sec, mpc_t x, mpc_t value);
  
  /* Function in getopts.c */
  void mps_parse_opts (mps_status * s, int argc, char *argv[]);
  mps_boolean mps_getopts (mps_opt ** opt, int *argc_ptr, char ***argv_ptr,
                           const char *opt_format);



/* functions in sort.c */
  void mps_fsort (mps_status * s);
  void mps_dsort (mps_status * s);
  void mps_msort (mps_status * s);

/* functions in solve.c */
  void mps_update (mps_status * s);
  void mps_fsrad (mps_status * s, mps_cluster * cluster, cplx_t sc, double *sr);
  void mps_dsrad (mps_status * s, mps_cluster * cluster, cdpe_t sc, rdpe_t sr);
  void mps_msrad (mps_status * s, mps_cluster * cluster, mpc_t sc, rdpe_t sr);

  mps_boolean mps_check_stop (mps_status * s);
  void mps_fsolve (mps_status * s, mps_boolean * d_after_f);
  void mps_dsolve (mps_status * s, mps_boolean d_after_f);
  void mps_msolve (mps_status * s);
  void mps_fpolzer (mps_status * s, int *it, mps_boolean * excep);
  void mps_dpolzer (mps_status * s, int *it, mps_boolean * excep);
  void mps_mpolzer (mps_status * s, int *it, mps_boolean * excep);

  /* Functions in modify.c */
  void mps_fmodify (mps_status * s, mps_boolean track_new_cluster);
  void mps_dmodify (mps_status * s, mps_boolean track_new_cluster);
  void mps_mmodify (mps_status * s, mps_boolean track_new_cluster);


/* functions in starting.c */
  double mps_maximize_distance (mps_status * s, double last_sigma,
                                mps_cluster_item * cluster, int n);
  void mps_fstart (mps_status * s, int n, mps_cluster_item * cluster, double clust_rad,
                   double g, rdpe_t eps_out, double fap[]);
  void mps_dstart (mps_status * s, int n, mps_cluster_item * cluster, rdpe_t clust_rad,
                   rdpe_t g, rdpe_t eps_out, rdpe_t dap[]);
  void mps_mstart (mps_status * s, int n, mps_cluster_item * cluster, rdpe_t clust_rad,
                   rdpe_t g, rdpe_t dap[], mpc_t gg);
  void mps_frestart (mps_status * s);
  void mps_drestart (mps_status * s);
  void mps_mrestart (mps_status * s);
  void mps_fshift (mps_status * s, int m, mps_cluster_item * cluster, double clust_rad,
                   cplx_t g, rdpe_t eps);
  void mps_dshift (mps_status * s, int m, mps_cluster_item * cluster, rdpe_t clust_rad,
                   cdpe_t g, rdpe_t eps);
  void mps_mshift (mps_status * s, int m, mps_cluster_item * cluster, rdpe_t clust_rad,
                   mpc_t g);

/* functions in stio.c */
  void mps_readroots (mps_status * s);
  void mps_countroots (mps_status * s);
  void mps_outroot (mps_status * s, int i, int num);
  void mps_output (mps_status * s);
  void mps_copy_roots (mps_status * s);
  void mps_dump (mps_status * s);
  void mps_dump_cluster_structure (mps_status * s, FILE * outstr);
  mps_boolean mps_is_a_tty (FILE * stream);
  void mps_warn (mps_status * st, char *s);
  void mps_error (mps_status * st, int args, ...);



/* functions in test.c */
  mps_boolean mps_inclusion (mps_status * s);

/* functions in cluster.c */
  void mps_cluster_detach (mps_status * s, mps_cluster * cluster);

/* functions in touch.c */
  mps_boolean mps_ftouchnwt (mps_status * s, double * frad, int n, int i, int j);
  mps_boolean mps_dtouchnwt (mps_status * s, int n, int i, int j);
  mps_boolean mps_mtouchnwt (mps_status * s, int n, int i, int j);
  mps_boolean mps_ftouchreal (mps_status * s, int n, int i);
  mps_boolean mps_dtouchreal (mps_status * s, int n, int i);
  mps_boolean mps_mtouchreal (mps_status * s, int n, int i);
  mps_boolean mps_ftouchimag (mps_status * s, int n, int i);
  mps_boolean mps_dtouchimag (mps_status * s, int n, int i);
  mps_boolean mps_mtouchimag (mps_status * s, int n, int i);
  mps_boolean mps_ftouchunit (mps_status * s, int n, int i);
  mps_boolean mps_dtouchunit (mps_status * s, int n, int i);
  mps_boolean mps_mtouchunit (mps_status * s, int n, int i);

/* functions in user.c */
  void mps_fnewton_usr (mps_status * st, cplx_t x, double *rad, cplx_t corr,
                        mps_boolean * again);
  void mps_dnewton_usr (mps_status * st, cdpe_t x, rdpe_t rad, cdpe_t corr,
                        mps_boolean * again);
  void mps_mnewton_usr (mps_status * st, mpc_t x, rdpe_t rad, mpc_t corr,
                        mps_boolean * again);

/* Routines of Input/Output in stio.c */
  void mps_skip_comments (FILE * input_stream);

  mps_input_option
  mps_parse_option_line (mps_status * s, char *line, size_t length);

  void
  mps_parse_stream (mps_status * s, FILE * input_stream);

  /* Functions in horner.c */
  void mps_fhorner (mps_status * s, mps_monomial_poly * p, cplx_t x, cplx_t value);
  void mps_fhorner_with_error (mps_status * s, mps_monomial_poly * p, cplx_t x, cplx_t value, double * relative_error);
  void mps_dhorner (mps_status * s, mps_monomial_poly * p, cdpe_t x, cdpe_t value);
  void mps_dhorner_with_error (mps_status * s, mps_monomial_poly * p, cdpe_t x, cdpe_t value, rdpe_t relative_error);
  void mps_mhorner (mps_status * s, mps_monomial_poly * p, mpc_t x, mpc_t value);
  void mps_mhorner_with_error (mps_status * s, mps_monomial_poly * p, mpc_t x, mpc_t value, rdpe_t relative_error, long int wp);
  void mps_mhorner_with_error2 (mps_status * s, mps_monomial_poly * p, mpc_t x, mpc_t value, rdpe_t relative_error, long int wp);

/*
 * End of extern "C" {
 *   ...
 * }
 */
#ifdef __cplusplus
}
#endif

#endif                          /* ndef MPS_CORE_H */
