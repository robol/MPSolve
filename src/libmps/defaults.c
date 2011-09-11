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

#include <float.h>
#include <mps/core.h>
#include <mps/threading.h>
#include <string.h>

void
mps_set_default_values (mps_status * s)
{

  /* flags */
  s->skip_float = false;        /* set to true to skip float phase     */
  s->resume = false;            /* resume from pre-computed roots      */
  s->chkrad = false;            /* check radii after completion        */

  /* I/O flags */
  s->DOLOG = false;             /* if nonzero enables the logging      */
  s->DOWARN = true;             /* if nonzero enables the warnings     */
  s->DOSORT = true;             /* if nonzero enables root sorting     */
  s->debug_level = 0;

  /* I/O streams */
  s->instr = NULL;              /* input stream                        */
  s->outstr = NULL;             /* output stream                       */
  s->logstr = NULL;             /* log stream                          */
  s->rtstr = NULL;              /* root stream                         */

  /* constants/parameters */
  s->max_pack = 1000;           /* number of max packets of iterations */
  s->max_it = 10;               /* number of max iterations per packet */
  s->max_newt_it = 15;          /* number of max newton iterations for */

  /* Set number of threads to 1.5 * number_of_cores, if this is
   * computable. Set it to 12 otherwise.                     */
  s->n_threads = (int) 1.5 *mps_thread_get_core_number (s);
  if (!s->n_threads)
    s->n_threads = 12;


  s->mpwp_max = 1000000000;     /* maximum allowed bits for mp         */
  /*   numbers: used in hi-prec. shifts  */
  strncpy (s->goal, "iannc", 5);        /* stores the goal, by default "iannc" */
  s->prec_out = 2 * DBL_DIG;    /* digits of required output precision */

  /* polynomial data - shared variables */
  s->prec_in = -1;              /* number of digits of input precision */
  /*   override input file if != -1      */
  s->data_type = NULL;          /* stores the input data type          */
  s->spar = NULL;               /* sparsity structure of polynomial    */
  s->fpr = NULL;                /* standard real coefficients          */
  s->fpc = NULL;                /* standard complex coefficients       */
  s->dpr = NULL;                /* dpe real coefficients               */
  s->dpc = NULL;                /* dpe complex coefficients            */
  s->mip_r = NULL;              /* real part of integer input coefs    */
  s->mip_i = NULL;              /* imag. part of integer input coefs   */
  s->mqp_r = NULL;              /* real part of rational input coeff.  */
  s->mqp_i = NULL;              /* imag. part of rational input coefs  */
  s->mfpr = NULL;               /* multiprecision real coefficients    */
  s->mfpc = NULL;               /* multiprecision complex coefficients */
  s->zero_roots = 0;

  /* soution related variables */
  s->lastphase = no_phase;      /* store last computed phase           */
  s->order = NULL;              /* output index order: ord[0],..,ord[n] */
  s->froot = NULL;              /* root approxs. as standard complex   */
  s->droot = NULL;              /* root approx. as complex dpe         */
  s->mroot = NULL;              /* root approximations as complex mp   */
  s->frad = NULL;               /* radii of inclusion disks as real    */
  s->drad = NULL;               /* radii of inclusion disks as rdpe_t  */

  s->clust = NULL;              /* indices of components of clusters   */
  s->punt = NULL;               /* beginning of each cluster           */
  s->rootwp = NULL;             /* working precision used for each root */
  s->mfppc = NULL;              /* multiprecision complex coeffs of p' */
  s->fap = NULL;                /* moduli of the coefficients as double */
  s->dap = NULL;                /* moduli of the coefficients as dpe   */
  s->again = NULL;              /* flag vector: true where more        */
  s->fppc = NULL;               /* standard complex coefficients       */
  s->fppc1 = NULL;              /* standard complex coefficients       */
  s->dpc1 = NULL;               /* dpe complex coefficients            */
  s->dpc2 = NULL;               /* dpe complex coefficients            */
  s->mfpc1 = NULL;              /* temp multiprec. complex coeff. of p' */
  s->mfpc2 = NULL;              /* temp multiprec. complex coeff.      */
  s->mfppc1 = NULL;             /* temp multiprec. complex coeff.      */

  s->spar1 = NULL;              /* temp sparsity structure of poly     */
  s->spar2 = NULL;              /* temp sparsity structure of poly     */
  s->oldpunt = NULL;            /* stores the previous value of punt   */
  s->fap1 = NULL;               /* moduli of the coefficients as double */
  s->fap2 = NULL;               /* temp. log of the coeffs as double   */
  s->dap1 = NULL;               /* temp moduli of the coeffs as dpe    */
  s->dap2 = NULL;               /* temp. log of the coeffs as double   */

  s->h = NULL;                  /* needed for convex hull computations */
  s->again_old = NULL;          /* temp flag vector: true where more   */
  /* iterations must be performed        */
  s->clust_aux = NULL;          /* auxiliary vectors ...               */
  s->punt_aux = NULL;
  s->punt_out = NULL;
  s->clust_out = NULL;

  /* Don't use user define functions in the default case */
  s->fnewton_usr = NULL;
  s->dnewton_usr = NULL;
  s->mnewton_usr = NULL;
  s->check_data_usr = NULL;
  s->fstart_usr = NULL;
  s->dstart_usr = NULL;

  s->random_seed = 0;

  /* Default algorithm */
  mps_status_select_algorithm (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
}
