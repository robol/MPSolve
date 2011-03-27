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
 * this file is needed to access all the MPSolve routines. 
 *
 * @brief Header file for libmps
 */

#ifndef __MPS_H__
#define __MPS_H__

/* standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <limits.h>

/* gmp library */
#include <gmp.h>

/* local include files */
#include "tools.h"
#include "mt.h"
#include "gmptools.h"
#include "mpc.h"
#include "link.h"
#include "mptemp.h"

/* constants */

/**
 * @brief Type representing the computation phase 
 * of the algorithm we are in 
 * now. It can assume the values:
 * - <code>no_phase</code>;
 * - <code>float_phase</code>;
 * - <code>dpe_phase</code>;
 * - <code>mp_phase</code>;
 */
typedef enum {no_phase, float_phase, dpe_phase, mp_phase} mps_phase;

/**
 * @brief this struct holds the state of the mps computation
 */
typedef struct {

   boolean resume;		/* to complete                         */
   boolean chkrad;		/* check radii after completion        */

   int newtis;
   int newtis_old;

  /*
   * INPUT / OUTPUT STREAMS
   */

  /**
   * @brief <code>true</code> if log are needed. They will
   * be written to <code>logstr</code>
   *
   * @see logstr
   */
   boolean DOLOG;

  /**
   * @brief <code>true</code> if warning are needed. 
   */
   boolean DOWARN;

  /**
   * @brief <code>true</code> if root sorting is desired. It will
   * be performed with routines in <code>mps_sort.c</code>.
   */
   boolean DOSORT;

  /**
   * @brief Default input stream.
   */
   FILE * instr;
  
  /**
   * @brief Default output stream.
   */
   FILE * outstr;

  /**
   * @brief Default log stream
   */
   FILE * logstr;
  
   FILE * rtstr;
  
  /*
   * CONSTANT, PARAMETERS
   */
  
  /**
   * @brief number of max packets of iterations 
   */
   int max_pack;		

  /**
   * @brief number of max iterations per packet 
   */
   int max_it;		
  
  /**
   * @brief Number of max newton iterations for gravity center
   * computations.
   */
   int max_newt_it;

  /**
   * @brief Maximum allowed number of bits for mp numbers: used in 
   * high precision shift.
   */
   long int mpwp_max;	

  /**
   * @brief stores the goal of the computation
   *
   * <code>goal</code> is an array of 5 chars with this meaning:
   * - <code>goal[0]</code> can assume the following values:
   *   - <code>a</code> that means \b approximate;
   *   - <code>c</code> that means \b count the roots;
   *   - <code>i</code> that means \b isolate the roots;
   * - <code>goal[1]</code>represents the search set for
   *   the roots. It can be:
   *   - <code>a</code> that means \f$ \mathbb{C} \f$;
   *   - <code>r</code> that means \f$ \{ z \in \mathbb{C} \ | \ \mathrm{Re}{z} > 0 \} \f$;
   *   - <code>l</code> that means \f$ \{ z \in \mathbb{C} \ | \ \mathrm{Re}{z} < 0 \} \f$;
   *   - <code>u</code> that means \f$ \{ z \in \mathbb{C} \ | \ \mathrm{Im}{z} > 0 \} \f$;
   *   - <code>d</code> that means \f$ \{ z \in \mathbb{C} \ | \ \mathrm{Im}{z} < 0 \} \f$;
   *   - <code>i</code> that means \f$ \{ z \in \mathbb{C} \ | \ \lVert z \rVert \leq 1 \} \f$;
   *   - <code>i</code> that means \f$ \{ z \in \mathbb{C} \ | \ \lVert z \rVert \geq 1 \} \f$;
   * - <code>goal[2]</code> controls if multiplicity check is enabled or not. Its value can be:
   *   - <code>m</code> if it is enabled;
   *   - <code>n</code> if that's not the case;
   * - <code>goal[3]</code> sets what properties of the roots have to be detected
   *   during the computation:
   *   - <code>n</code> means no one, and this is the default;
   *   - <code>r</code> means detect if a root is real;
   *   - <code>i</code> means detect if a root is imaginary;
   *   - <code>b</code> means both the above;
   * - <code>goal[4]</code> determines output format for the results:
   *   - <code>b</code> means \b bare and it's the default value;
   *   - <code>g</code> means \b gnuplot format;
   *   - <code>c</code> means \b compact;
   *   - <code>v</code> means \b verbose;
   *   - <code>f</code> means \b full;
   */
   char goal[5];
  
  /**
   * @brief digits of required output precision
   *
   */
   long int prec_out;
  
  /**
   * @brief boolean value that determine if we should
   * use a random seed for startingd points
   */
   boolean random_seed;
  
  /*
   * POLYNOMIAL DATA: SHARED VARIABLES
   */
  
  /**
   * @brief degree of zero-deflated polynomial.
   */
   int n;
  
  /**
   * @brief input degree and allocation size.
   */
   int deg;
  
  /**
   * @brief stores the input data type
   *
   * The value of <code>data_type[0]</code> can be:
   * - <code>'s'</code> if the input is a sparse polynomial;
   * - <code>'u'</code> if the input is a user polynomial;
   * - <code>'d'</code> if the input is a dense polybomial;
   *
   * while the value of <code>data_type[1]</code> can be:
   * - <code>'r'</code> if the coefficients are real;
   * - <code>'c'</code> if the coefficents are complex;
   *
   * and finally <code>data_type[2]</code> can assume the following values:
   * - <code>'i'</code> means integer coefficients;
   * - <code>'q'</code> means rationa cofficients;
   * - <code>'b'</code> means bigfloat coefficents;
   * - <code>'f'</code> means floating point coefficients;
   */
   char * data_type;
  
  /**
   * @brief Number of digits of input precision in its binary
   * representation. 
   */
   long int prec_in;
  
  /**
   * @brief This array contains the structure of the sparse
   * polynomial. 
   *
   * <code>spar[i]</code> is <code>true</code> if and only if
   * the i-th coefficients of the polynomial is a non-zero 
   * coefficients
   */
   boolean * spar;		/* sparsity structure of the poly.     */
  
   double * fpr;		/* standard real coefficients          */
   cplx_t * fpc;		/* standard complex coefficients       */
   rdpe_t * dpr;		/* dpe real coefficients               */
   cdpe_t * dpc;		/* dpe complex coefficients            */
   mpz_t * mip_r;		/* real part of integer input coeffs   */
   mpz_t * mip_i;		/* imag. part of integer input coeffs  */
   mpq_t * mqp_r;		/* real part of rational input coeffs  */
   mpq_t * mqp_i;		/* imag. part of rational input coeffs */
   mpf_t * mfpr;		/* multiprecision real coefficients    */
   mpc_t * mfpc;		/* multiprecision complex coefficients */
  
  /* soution related variables */
   mps_phase lastphase;		/* store last computed phase           */
  
  /**
   * @brief Vector containing count of in, out and uncertaing roots.
   */
   int count[3];
  
  /**
   * @brief Number of zero roots.
   */
   int zero_roots;		/* number of roots = 0                 */
  
  /**
   * @brief Status of each approximation
   * 
   * <code>status</code> is an array of char arrays
   * composed of three elements. For every <code>i</code>:
   * - <code>status[i][0]</code> can assume the values:
   *   - <code>m</code>: multiple root;
   *   - <code>i</code>: isolated root;
   *   - <code>a</code>: approximated single root (relative error
   *     is less then \f$ 10^{-d_{out}} \f$);
   *   - <code>o</code>: approximated cluster of roots (relative
   *     error is less than \f$ 10^{-d_{out}} \f$);
   *   - <code>c</code>: cluster of roots not yet approximated (relative
   *     error is greater than \f$ 10^{-d_{out}} \f$);
   *   - <code>f</code>: TODO: Determine what this state means;
   *   - <code>x</code>: this root cannot be represented as a <code>double</code>, i.e. it
   *     is <code>< DBL_MIN</code> or <code>> DBL_MAX</code>;
   * - <code>status[i][1]</code> can assume the values:
   *   - <code>R</code>: real root;
   *   - <code>r</code>: non real root;
   *   - <code>I</code>: imaginary root;
   *   - <code>i</code>: non imaginary root;
   *   - <code>w</code>: uncertain real/imaginary root;
   *   - <code>z</code>: non real and non imaginary root;
   * - <code>status[i][2]</code> can assume the values:
   *   - <code>i</code>: root in \f$ \mathcal{S} \f$;
   *   - <code>o</code>: root out of \f$ \mathcal{S} \f$;
   *   - <code>u</code>: root uncertain;
   */
   char (* status)[3];	/* status of each approximation        */
   int * order;		/* output index order: ord[0],..,ord[n]*/
   cplx_t * froot;		/* root approx. as float complex numbers */
   cdpe_t * droot;		/* root approx. as complex dpe numbers */
   mpc_t * mroot;		/* root approx. as complex mp numbers  */
   double * frad;		/* radii of incl. disks as real nums   */
   rdpe_t * drad;		/* radii of incl. disks as rdpe_t nums */
  
  /* lifetime global variables */
   boolean skip_float;	/* skip float phase?                   */
   rdpe_t eps_in;		/* input precision                     */
   rdpe_t eps_out;		/* output precision                    */
   double lmax_coeff;	/* log of the max modulus of coeffs.   */
   long int mpwp;		/* # of bits of current working prec.  */
   rdpe_t mp_epsilon;	/* current mp epsilon                  */
   double sep;		/* Log of the lower bound to the 
                                   minimum distance of the roots       */
  
  /**
   * @brief number of active clusters
   */
   int nclust;		


  /**
   * @brief indices of cluster components
   *
   * <code>clust</code> is an integer array containing the indexes
   * of roots in every cluster. 
   *
   * The indexes of the <code>j</code>th cluster are
   * <code>clust[punt[j] : punt[j+1]]
   *
   * @see punt
   */
   int * clust;
  
  /**
   * @brief begginning of each cluster
   *
   * <code>punt</code> is a vector of <code>nclust</code> integers;
   *
   * For every <code>i</code> <code>punt[i]</code> is the index in
   * the integer vector <code>clust</code> corresponding to the first
   * index of a cluster, i.e. the jth cluster of roots is composed by
   * roots indexed on <code>clust[p[j] : p[j+1]]</code>.
   *
   * @see nclust
   * @see clust
   */
   int * punt;
  
  /**
   * @brief Array containing working precisions used for each root.
   */
   long int * rootwp;
  
  /**
   * @brief Multiprecision complex coefficients of \f$p'(x)\f$.
   */
   mpc_t * mfppc;		
  
  /**
   * @brief Array containing moduli of the coefficients as double numbers.
   */
   double * fap;
  
  /**
   * @brief Array containing moduli of the coefficients as dpe numbers.
   */
   rdpe_t * dap;
  
  /**
   * @brief Array that whose i-th component is set to <code>true</code> if
   * the i-th root needs more iterations.
   */
   boolean * again;		
  
  /**
   * @brief Array containing standard complex coefficients
   */
   cplx_t * fppc;
  
   cplx_t * fppc1;		/* standard complex coefficients       */
   cdpe_t * dpc1;		/* dpe complex coefficients            */
   cdpe_t * dpc2;		/* dpe complex coefficients            */
   mpc_t * mfpc1;		/* temp multiprecision complex coeffs  */
   mpc_t * mfpc2;		/* temp multiprecision complex coeffs  */
   mpc_t * mfppc1;		/* temp multiprec complex coeffs of p' */
  
   boolean * spar1;		/* temp sparsity structure of the poly */
   boolean * spar2;		/* temp sparsity structure of the poly */
   int * oldpunt;		/* stores the previous value of punt   */
   double * fap1;		/* moduli of the coeffs as double      */
   double * fap2;		/* log of the moduli of the coeff. as double */
   rdpe_t * dap1;		/* moduli of the coeff.as dpe numbers  */
   rdpe_t * dap2;		/* moduli of the log of the coeff. as double */
  
  /**
   * @brief Vector needed for convex hull computation. 
   *
   * It is <code>true</code> in position \f$j\f$ if 
   * the point \f$(j, log(x_j))\f$ is a vertex of the convex
   * hull computed by <code>fconvex()</code> and the other functions in
   * <code>mps_cnvx.c</code>
   */
   boolean * h;
  
   boolean * again_old;	/* temp flag vector: true in the  components
                                   where iterations must be performed  */
  
   int * clust_aux;		/* auxiliary vector                    */
   int * punt_aux;		/* auxiliary vector                    */
   int * punt_out;		/* auxiliary vector                    */
   int * clust_out;		/* auxiliary vector                    */

} mps_status; /* End of typedef struct { ... */


/* FUNCTIONS */

/* functions in mps_defaults.c */
void mps_set_default_values(mps_status* s);

/* functions in mps_aber.c */
void mps_faberth(mps_status* s, int j, cplx_t abcorr);
void mps_daberth(mps_status* s, int j, cdpe_t abcorr);
void mps_maberth(mps_status* s, int j, mpc_t abcorr);
void mps_faberth_s(mps_status* s, int j, int jc, cplx_t abcorr);
void mps_daberth_s(mps_status* s, int j, int jc, cdpe_t abcorr);
void mps_maberth_s(mps_status* s, int j, int jc, mpc_t abcorr);
void mps_mnewtis(mps_status* s);

/* functions in mps_clus.c */
void mps_fcluster(mps_status* s, int nf);
void mps_dcluster(mps_status* s, int nf);
void mps_mcluster(mps_status* s, int nf);

/* functions in mps_cnvx.c */
void mps_fconvex(int n, double a[]);

/* functions in mps_data.c */
void mps_mp_set_prec(mps_status* s, long int prec);
void mps_allocate_data(mps_status* s);
void mps_prepare_data(mps_status* s, long int prec);
void mps_restore_data(mps_status* s);
void mps_free_data(mps_status* s);

/* functions in mps_impr.c */
void mps_improve(void);

/* functions in mps_main.c */
void mps_mpsolve(mps_status* s);
void mps_setup(mps_status* s);
void mps_check_data(mps_status* s, char * which_case);
void mps_compute_sep(mps_status* s);

/* functions in mps_newt.c */
void mps_fnewton(int n, cplx_t z, double * radius, cplx_t corr, cplx_t fpc[],
		 double fap[], boolean *  cont);
void mps_dnewton(int n, cdpe_t z, rdpe_t radius, cdpe_t corr,
		 cdpe_t dpc[], rdpe_t dap[], boolean * cont);
void mps_mnewton(int n, mpc_t z, rdpe_t radius, mpc_t corr, mpc_t mfpc[],
		 mpc_t mfppc[], rdpe_t dap[], boolean * spar, boolean * cont);
void mps_parhorner(int n, mpc_t x, mpc_t p[], boolean b[], mpc_t s);
void mps_aparhorner(int n, rdpe_t x, rdpe_t p[], boolean b[], rdpe_t s);

/**
 * This function parse command lines and sets global variables 
 * defined in mps.h in an appropriate way. 
 *
 * It is implemented in mps_opts.c
 *
 * @brief Parse options from command line.
 *
 * @param argc 
 *   Argoment counter as obtained from the main()
 *   function.
 * @param argv 
 *   Argoment values, i.e. array of char* as obtained
 *   in the main() function.
 */
void mps_parse_opts(mps_status* s, int argc, char *argv[]);

/* functions in mps_sort.c */
void mps_fsort(mps_status* s);
void mps_dsort(mps_status* s);
void mps_msort(mps_status* s);

/* functions in mps_solv.c */
void mps_update(void);
void mps_fsrad(int i, cplx_t sc, double *sr);
void mps_dsrad(int i, cdpe_t sc, rdpe_t sr);
void mps_msrad(int i, mpc_t sc, rdpe_t sr);
void mps_fmodify(void);
void mps_dmodify(void);
void mps_mmodify(void);
boolean mps_check_stop();
void mps_fsolve(boolean * d_after_f);
void mps_dsolve(boolean d_after_f);
void mps_msolve(void);
void mps_fpolzer(int * it, boolean * excep);
void mps_dpolzer(int * it, boolean * excep);
void mps_mpolzer(int * it, boolean * excep);

/* functions in mps_star.c */
void mps_fstart(int n, int i_clust, double clust_rad, double g, rdpe_t eps_out,
		double fap[]);
void mps_dstart(int n, int i_clust, rdpe_t clust_rad, rdpe_t g, rdpe_t eps_out,
	    rdpe_t dap[]);
void mps_mstart(int n, int i_clust, rdpe_t clust_rad, rdpe_t g, rdpe_t dap[]);
void mps_frestart(void);
void mps_drestart(void);
void mps_mrestart(void);
void mps_fshift(int m, int i_clust, double clust_rad, cplx_t g, rdpe_t eps);
void mps_dshift(int m, int i_clust, rdpe_t clust_rad, cdpe_t g, rdpe_t eps);
void mps_mshift(int m, int i_clust, rdpe_t clust_rad, mpc_t g);

/* functions in mps_stio.c */
void mps_readroots(void);
void mps_countroots(void);
void mps_outroot(int i);
void mps_output(void);
void mps_copy_roots(void);
void mps_dump(FILE * dmpstr);
void mps_warn(char * s);
void mps_error(int args, ...);

/* functions in mps_test.c */
boolean  mps_inclusion(void);

/* functions in mps_touch.c */
boolean mps_ftouchnwt(mps_status* s, int n, int i, int j);
boolean mps_dtouchnwt(mps_status* s, int n, int i, int j);
boolean mps_mtouchnwt(mps_status* s, int n, int i, int j);
boolean mps_ftouchreal(mps_status* s, int n, int i);
boolean mps_dtouchreal(mps_status* s, int n, int i);
boolean mps_mtouchreal(mps_status* s, int n, int i);
boolean mps_ftouchimag(mps_status* s, int n, int i);
boolean mps_dtouchimag(mps_status* s, int n, int i);
boolean mps_mtouchimag(mps_status* s, int n, int i);
boolean mps_ftouchunit(mps_status* s, int n, int i);
boolean mps_dtouchunit(mps_status* s, int n, int i);
boolean mps_mtouchunit(mps_status* s, int n, int i);

/* functions in mps_usr.c */
void mps_fnewton_usr(cplx_t x, double * rad, cplx_t corr, boolean * again);
void mps_dnewton_usr(cdpe_t x, rdpe_t rad, cdpe_t corr, boolean * again);
void mps_mnewton_usr(mpc_t x, rdpe_t rad, mpc_t corr, boolean * again);

#endif /* ndef __MPS_H__ */



