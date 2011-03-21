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

#include "mps.h"

/* global data declaration and initialization                          */

/* flags */
boolean skip_float = false;	/* set to true to skip float phase     */
boolean resume = false; 	/* resume from pre-computed roots      */
boolean chkrad = false;		/* check radii after completion        */
int newtis;                     /* check for newton-isolation          */
int newtis_old;                 /* store old value of newtis           */

/* I/O flags */
boolean DOLOG = false;		/* if nonzero enables the logging      */
boolean DOWARN = true;		/* if nonzero enables the warnings     */
boolean DOSORT = true;		/* if nonzero enables root sorting     */

/* I/O streams */
FILE *instr = NULL;		/* input stream                        */
FILE *outstr = NULL;		/* output stream                       */
FILE *logstr = NULL;		/* log stream                          */
FILE *rtstr = NULL;		/* root stream                         */

/* constants/parameters */
int max_pack = 1000;		/* number of max packets of iterations */
int max_it = 10;		/* number of max iterations per packet */
int max_newt_it = 15;		/* number of max newton iterations for */
				/*   gravity center computations       */
long int mpwp_max = 1000000000;	/* maximum allowed bits for mp         */ 
				/*   numbers: used in hi-prec. shifts  */
char goal[5] = "iannc";		/* stores the goal, by default "iannc" */
long int prec_out = 2*DBL_DIG;  /* digits of required output precision */

/* polynomial data - shared variables */
int n = 0;			/* degree of zero-deflated polynomial  */
int deg = 0;			/* input degree and allocation size    */
long int prec_in = -1;		/* number of digits of input precision */
				/*   override input file if != -1      */
boolean random_seed;            /* true if random seed should be used  */
char *data_type = NULL;		/* stores the input data type          */
boolean *spar = NULL;		/* sparsity structure of polynomial    */
double *fpr = NULL;		/* standard real coefficients          */
cplx_t *fpc = NULL;		/* standard complex coefficients       */
rdpe_t *dpr = NULL;		/* dpe real coefficients               */
cdpe_t *dpc = NULL;		/* dpe complex coefficients            */
mpz_t *mip_r = NULL;		/* real part of integer input coefs    */
mpz_t *mip_i = NULL;		/* imag. part of integer input coefs   */
mpq_t *mqp_r = NULL;		/* real part of rational input coeff.  */
mpq_t *mqp_i = NULL;		/* imag. part of rational input coefs  */
mpf_t *mfpr = NULL;		/* multiprecision real coefficients    */
mpc_t *mfpc = NULL;		/* multiprecision complex coefficients */

/* soution related variables */
phase lastphase = no_phase;	/* store last computed phase           */
int count[3];			/* count roots: [in, out, uncertain]   */
int zero_roots;			/* number of roots = 0                 */
char (*status)[3];		/* status of each approximation        */
int *order = NULL;		/* output index order: ord[0],..,ord[n]*/
cplx_t *froot = NULL;		/* root approxs. as standard complex   */
cdpe_t *droot = NULL;		/* root approx. as complex dpe         */
mpc_t *mroot = NULL;		/* root approximations as complex mp   */
double *frad = NULL;		/* radii of inclusion disks as real    */
rdpe_t *drad = NULL;		/* radii of inclusion disks as rdpe_t  */

/* lifetime global variables */
rdpe_t eps_in;			/* input precision                     */
rdpe_t eps_out;			/* output precision                    */
double lmax_coeff;		/* log of the max modulus of coeffs.   */
long int mpwp;			/* bits of current working precision   */
rdpe_t mp_epsilon;		/* current mp epsilon                  */
double sep;			/* Log of the lower bound to the       */
				/* miniumum distance of the roots      */
int nclust;			/* number of active clusters           */
int *clust = NULL;		/* indices of components of clusters   */
int *punt = NULL;		/* beginning of each cluster           */
long int *rootwp = NULL;	/* working precision used for each root*/
mpc_t *mfppc = NULL;		/* multiprecision complex coeffs of p' */
double *fap = NULL;		/* moduli of the coefficients as double*/
rdpe_t *dap = NULL;		/* moduli of the coefficients as dpe   */
boolean *again = NULL;		/* flag vector: true where more        */
				/* iterations must be performed        */
cplx_t *fppc = NULL;		/* standard complex coefficients       */
cplx_t *fppc1 = NULL;		/* standard complex coefficients       */
cdpe_t *dpc1 = NULL;		/* dpe complex coefficients            */
cdpe_t *dpc2 = NULL;		/* dpe complex coefficients            */
mpc_t *mfpc1 = NULL;		/* temp multiprec. complex coeff. of p'*/
mpc_t *mfpc2 = NULL;		/* temp multiprec. complex coeff.      */
mpc_t *mfppc1 = NULL;		/* temp multiprec. complex coeff.      */

boolean *spar1 = NULL;		/* temp sparsity structure of poly     */
boolean *spar2 = NULL;		/* temp sparsity structure of poly     */
int *oldpunt = NULL;		/* stores the previous value of punt   */
double *fap1 = NULL;		/* moduli of the coefficients as double*/
double *fap2 = NULL;		/* temp. log of the coeffs as double   */
rdpe_t *dap1 = NULL;		/* temp moduli of the coeffs as dpe    */
rdpe_t *dap2 = NULL;		/* temp. log of the coeffs as double   */

boolean *h = NULL;		/* needed for convex hull computations */
boolean *again_old = NULL;	/* temp flag vector: true where more   */
				/* iterations must be performed        */
int *clust_aux = NULL;          /* auxiliary vectors ...               */
int *punt_aux = NULL;
int *punt_out = NULL;
int *clust_out = NULL;
