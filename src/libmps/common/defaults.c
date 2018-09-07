/*
 * This file is part of MPSolve 3.1.6
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@sns.it>
 */


#include <stdio.h>
#include <float.h>
#include <mps/mps.h>
#include <string.h>

MPS_PRIVATE void
mps_set_default_values (mps_context * s)
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
  s->instr = stdin;              /* input stream                        */
  s->outstr = stdout;             /* output stream                       */
  s->logstr = stderr;             /* log stream                          */
  s->rtstr = NULL;              /* root stream                         */

  /* constants/parameters */
  s->max_pack = 100000;           /* number of max packets of iterations */
  s->max_it = 20;                /* number of max iterations per packet */
  s->max_newt_it = 15;           /* number of max newton iterations for */
  s->jacobi_iterations = false;

  /* Set number of threads to 1.5 * number_of_cores, if this is
   * computable. Set it to 12 otherwise.                     */
  s->n_threads = (int)1.5 * mps_thread_get_core_number (s);
  if (!s->n_threads)
    s->n_threads = 12;

  s->clusterization = NULL;

  s->mpwp_max = 100000000;     /* maximum allowed bits for mp         */

  s->zero_roots = 0;

  /* soution related variables */
  s->lastphase = no_phase;      /* store last computed phase           */
  s->order = NULL;              /* output index order: ord[0],..,ord[n] */
  s->fppc1 = NULL;              /* standard complex coefficients       */
  s->dpc1 = NULL;               /* dpe complex coefficients            */
  s->dpc2 = NULL;               /* dpe complex coefficients            */
  s->mfpc1 = NULL;              /* temp multiprec. complex coeff. of p' */
  s->mfppc1 = NULL;             /* temp multiprec. complex coeff.      */

  s->spar1 = NULL;              /* temp sparsity structure of poly     */
  s->oldpunt = NULL;            /* stores the previous value of punt   */
  s->fap1 = NULL;               /* moduli of the coefficients as double */
  s->fap2 = NULL;               /* temp. log of the coeffs as double   */
  s->dap1 = NULL;               /* temp moduli of the coeffs as dpe    */

  s->again_old = NULL;          /* temp flag vector: true where more   */

  s->random_seed = 0;
  s->newtis = 0;
  s->last_sigma = 0.1;

  s->secular_equation = NULL;
  s->active_poly = NULL;

  /* Input */
  s->input_config->starting_phase = no_phase;

  /* Output */
  s->output_config->format = MPS_OUTPUT_FORMAT_COMPACT;
  s->output_config->prec = 0.8 * DBL_DIG * LOG2_10;
  s->output_config->goal = MPS_OUTPUT_GOAL_ISOLATE;
  s->output_config->multiplicity = false;
  s->output_config->root_properties = MPS_OUTPUT_PROPERTY_NONE;
  s->output_config->search_set = MPS_SEARCH_SET_COMPLEX_PLANE;

  s->data_prec_max.value = 53;

  s->mpwp = DBL_DIG * LOG2_10;

  /* Setup for the default algorithm */
  s->mpsolve_ptr = MPS_MPSOLVE_PTR (mps_standard_mpsolve);
  s->algorithm = MPS_ALGORITHM_STANDARD_MPSOLVE;
  s->starting_strategy = MPS_STARTING_STRATEGY_DEFAULT;

  /* Allocate the thread_pool used in computations. */
  s->pool = mps_thread_pool_new (s, 0);

  /* Callbacks for async version */
  s->callback = NULL;
  s->user_data = NULL;

  /* Error handling */
  s->error_state = false;
  s->last_error = NULL;

  s->over_max = false;

  s->bmpc = NULL;

  /* Set a sensible gnuplot_format by default so the user won't see (null) in the
   * output. */
  s->gnuplot_format = "points";

  s->self_thread_pool = NULL;
  s->avoid_multiprecision = false;
  s->crude_approximation_mode = false;

  /* Get a standard regeneration driver to use for MPSolve. Note that this functions
   * does not actually allocate anything, but simply provides a pointer to the internal
   * MPSolve implementation of the standard regeneration_driver. This instance must
   * not be freed, since it is statically declared inside 
   * secsolve/standard-regeneration-driver.c */   
  s->regeneration_driver = mps_regeneration_driver_new_standard (s);
}
