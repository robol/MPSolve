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
typedef enum {no_phase, float_phase, dpe_phase, mp_phase} phase;

extern boolean resume;		/* to complete                         */
extern boolean chkrad;		/* check radii after completion        */

extern int newtis;
extern int newtis_old;

/*
 * INPUT / OUTPUT STREAMS
 */

/**
 * @brief <code>true</code> if log are needed. The will
 * be written to <code>logstr</code>
 *
 * @see logstr
 */
extern boolean DOLOG;

/**
 * @brief <code>true</code> if warning are needed. 
 */
extern boolean DOWARN;

/**
 * @brief <code>true</code> if root sorting is desired. It will
 * be performed with routines in <code>mps_sort.c</code>.
 */
extern boolean DOSORT;

/**
 * @brief Default input stream.
 */
extern FILE * instr;

/**
 * @brief Default output stream.
 */
extern FILE * outstr;

/**
 * @brief Default log stream
 */
extern FILE * logstr;

extern FILE * rtstr;

/*
 * CONSTANT, PARAMETERS
 */

/**
 * @brief number of max packets of iterations 
 */
extern int max_pack;		

/**
 * @brief number of max iterations per packet 
 */
extern int max_it;		

/**
 * @brief Number of max newton iterations for gravity center
 * computations.
 */
extern int max_newt_it;

/**
 * @brief Maximum allowed number of bits for mp numbers: used in 
 * high precision shift.
 */
extern long int mpwp_max;	

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
extern char goal[5];

/**
 * @brief digits of required output precision
 *
 */
extern long int prec_out;

 /*
 * POLYNOMIAL DATA: SHARED VARIABLES
 */

/**
 * @brief degree of zero-deflated polynomial.
 */
extern int n;

/**
 * @brief input degree and allocation size.
 */
extern int deg;

/**
 * @brief stores the input data type
 *
 * The value of <code>data_type[0]</code> can be:
 * - <code>'s'</code> if the input is a sparse polynomial;
 * - <code>'u'</code> if the input is a user polynomial;
 * - <code>'d'</code> if the input is a dense polybomial;
 * while the value of <code>data_type[1]</code> can be:
 * - <code>'r'</code> if the coefficients are real;
 * - <code>'c'</code> if the coefficents are complex;
 * and finally <code>data_type[2]</code> can assume the following values:
 * - <code>'i'</code> means integer coefficients;
 * - <code>'q'</code> means rationa cofficients;
 * - <code>'b'</code> means bigfloat coefficents;
 * - <code>'f'</code> means floating point coefficients;
 */
extern char * data_type;

/**
 * @brief Number of digits of input precision in its binary
 * representation. 
 */
extern long int prec_in;

/**
 * @brief This array contains the structure of the sparse
 * polynomial. 
 *
 * <code>spar[i]</code> is <code>true</code> if and only if
 * the i-th coefficients of the polynomial is a non-zero 
 * coefficients
 */
extern boolean * spar;		/* sparsity structure of the poly.     */

extern double * fpr;		/* standard real coefficients          */
extern cplx_t * fpc;		/* standard complex coefficients       */
extern rdpe_t * dpr;		/* dpe real coefficients               */
extern cdpe_t * dpc;		/* dpe complex coefficients            */
extern mpz_t * mip_r;		/* real part of integer input coeffs   */
extern mpz_t * mip_i;		/* imag. part of integer input coeffs  */
extern mpq_t * mqp_r;		/* real part of rational input coeffs  */
extern mpq_t * mqp_i;		/* imag. part of rational input coeffs */
extern mpf_t * mfpr;		/* multiprecision real coefficients    */
extern mpc_t * mfpc;		/* multiprecision complex coefficients */

/* soution related variables */
extern phase lastphase;		/* store last computed phase           */
extern int count[3];		/* count roots: [in, out, uncertain]   */
extern int zero_roots;		/* number of roots = 0                 */

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
extern char (* status)[3];	/* status of each approximation        */
extern int * order;		/* output index order: ord[0],..,ord[n]*/
extern cplx_t * froot;		/* root approx. as float complex numbers */
extern cdpe_t * droot;		/* root approx. as complex dpe numbers */
extern mpc_t * mroot;		/* root approx. as complex mp numbers  */
extern double * frad;		/* radii of incl. disks as real nums   */
extern rdpe_t * drad;		/* radii of incl. disks as rdpe_t nums */

/* lifetime global variables */
extern boolean skip_float;	/* skip float phase?                   */
extern rdpe_t eps_in;		/* input precision                     */
extern rdpe_t eps_out;		/* output precision                    */
extern double lmax_coeff;	/* log of the max modulus of coeffs.   */
extern long int mpwp;		/* # of bits of current working prec.  */
extern rdpe_t mp_epsilon;	/* current mp epsilon                  */
extern double sep;		/* Log of the lower bound to the 
                                   minimum distance of the roots       */

/**
 * @brief number of active clusters
 */
extern int nclust;		


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
extern int * clust;

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
 *
 * @see nclust
 * @see clust
 */
extern int * punt;

extern long int * rootwp;	/* working prec. used for each root    */

extern mpc_t * mfppc;		/* multiprecision complex coeffs of p' */

extern double * fap;		/* moduli of the coeffs as double numbers */
extern rdpe_t * dap;		/* moduli of the coefficients as dpe numbers */
extern boolean * again;		/* flag vector: true in the components where
				 * iterations must be performed        */

extern cplx_t * fppc;		/* standard complex coefficients       */
extern cplx_t * fppc1;		/* standard complex coefficients       */
extern cdpe_t * dpc1;		/* dpe complex coefficients            */
extern cdpe_t * dpc2;		/* dpe complex coefficients            */
extern mpc_t * mfpc1;		/* temp multiprecision complex coeffs  */
extern mpc_t * mfpc2;		/* temp multiprecision complex coeffs  */
extern mpc_t * mfppc1;		/* temp multiprec complex coeffs of p' */

extern boolean * spar1;		/* temp sparsity structure of the poly */
extern boolean * spar2;		/* temp sparsity structure of the poly */
extern int * oldpunt;		/* stores the previous value of punt   */
extern double * fap1;		/* moduli of the coeffs as double      */
extern double * fap2;		/* log of the moduli of the coeff. as double */
extern rdpe_t * dap1;		/* moduli of the coeff.as dpe numbers  */
extern rdpe_t * dap2;		/* moduli of the log of the coeff. as double */

extern boolean * h;		/* needed for convex hull computations */
extern boolean * again_old;	/* temp flag vector: true in the  components
                                   where iterations must be performed  */

extern int * clust_aux;		/* auxiliary vector                    */
extern int * punt_aux;		/* auxiliary vector                    */
extern int * punt_out;		/* auxiliary vector                    */
extern int * clust_out;		/* auxiliary vector                    */


/* FUNCTIONS */


/* functions in mps_aber.c */
void faberth(int j, cplx_t abcorr);
void daberth(int j, cdpe_t abcorr);
void maberth(int j, mpc_t abcorr);
void faberth_s(int j, int jc, cplx_t abcorr);
void daberth_s(int j, int jc, cdpe_t abcorr);
void maberth_s(int j, int jc, mpc_t abcorr);
void mnewtis(void);

/* functions in mps_clus.c */
void fcluster(int nf);
void dcluster(int nf);
void mcluster(int nf);

/* functions in mps_cnvx.c */
void fconvex(int n, double a[]);

/* functions in mps_data.c */
void mp_set_prec(long int prec);
void allocate_data(void);
void prepare_data(long int prec);
void restore_data(void);
void free_data(void);

/* functions in mps_impr.c */
void improve(void);

/* functions in mps_main.c */
void mpsolve(void);
void setup(void);
void check_data(char * which_case);
void compute_sep(void);

/* functions in mps_newt.c */
void fnewton(int n, cplx_t z, double * radius, cplx_t corr, cplx_t fpc[],
	     double fap[], boolean *  cont);
void dnewton(int n, cdpe_t z, rdpe_t radius, cdpe_t corr,
	     cdpe_t dpc[], rdpe_t dap[], boolean * cont);
void mnewton(int n, mpc_t z, rdpe_t radius, mpc_t corr, mpc_t mfpc[],
	     mpc_t mfppc[], rdpe_t dap[], boolean * spar, boolean * cont);
void parhorner(int n, mpc_t x, mpc_t p[], boolean b[], mpc_t s);
void aparhorner(int n, rdpe_t x, rdpe_t p[], boolean b[], rdpe_t s);

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
void parse_opts(int argc, char *argv[]);

/* functions in mps_sort.c */
void fsort(void);
void dsort(void);
void msort(void);

/* functions in mps_solv.c */
void update(void);
void fsrad(int i, cplx_t sc, double *sr);
void dsrad(int i, cdpe_t sc, rdpe_t sr);
void msrad(int i, mpc_t sc, rdpe_t sr);
void fmodify(void);
void dmodify(void);
void mmodify(void);
boolean check_stop();
void fsolve(boolean * d_after_f);
void dsolve(boolean d_after_f);
void msolve(void);
void fpolzer(int * it, boolean * excep);
void dpolzer(int * it, boolean * excep);
void mpolzer(int * it, boolean * excep);

/* functions in mps_star.c */
void fstart(int n, int i_clust, double clust_rad, double g, rdpe_t eps_out,
	    double fap[]);
void dstart(int n, int i_clust, rdpe_t clust_rad, rdpe_t g, rdpe_t eps_out,
	    rdpe_t dap[]);
void mstart(int n, int i_clust, rdpe_t clust_rad, rdpe_t g, rdpe_t dap[]);
void frestart(void);
void drestart(void);
void mrestart(void);
void fshift(int m, int i_clust, double clust_rad, cplx_t g, rdpe_t eps);
void dshift(int m, int i_clust, rdpe_t clust_rad, cdpe_t g, rdpe_t eps);
void mshift(int m, int i_clust, rdpe_t clust_rad, mpc_t g);

/* functions in mps_stio.c */
void readroots(void);
void countroots(void);
void outroot(int i);
void output(void);
void copy_roots(void);
void dump(FILE * dmpstr);
void warn(char * s);
void error(int args, ...);

/* functions in mps_test.c */
boolean  inclusion(void);

/* functions in mps_touch.c */
boolean ftouchnwt(int n, int i, int j);
boolean dtouchnwt(int n, int i, int j);
boolean mtouchnwt(int n, int i, int j);
boolean ftouchreal(int n, int i);
boolean dtouchreal(int n, int i);
boolean mtouchreal(int n, int i);
boolean ftouchimag(int n, int i);
boolean dtouchimag(int n, int i);
boolean mtouchimag(int n, int i);
boolean ftouchunit(int n, int i);
boolean dtouchunit(int n, int i);
boolean mtouchunit(int n, int i);

/* functions in mps_usr.c */
void fnewton_usr(cplx_t x, double * rad, cplx_t corr, boolean * again);
void dnewton_usr(cdpe_t x, rdpe_t rad, cdpe_t corr, boolean * again);
void mnewton_usr(mpc_t x, rdpe_t rad, mpc_t corr, boolean * again);

#endif /* ndef __MPS_H__ */



