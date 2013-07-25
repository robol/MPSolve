/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#define _GNU_SOURCE

#include <mps/mps.h>
#include <string.h>

long int
mps_context_get_data_prec_max (mps_context * s)
{
  long int ret;
  MPS_LOCK (s->data_prec_max);
  ret = s->data_prec_max.value;
  MPS_UNLOCK (s->data_prec_max);
  return ret;
}

int 
mps_context_get_degree (mps_context * s)
{
  return s->n;
}

/**
 * @brief Select algorithm to use for computation.
 *
 * Valid values for this field are
 * - MPS_ALGORITHM_STANDARD_MPSOLVE for the standard MPSolve algorithm;
 * - MPS_ALGORITHM_SECULAR_MPSOLVE  for the standard MPSolve algorithm
 *   applied to secular equations;
 * - MPS_ALGORITHM_SECULAR_GA for the algorithm based on coefficient regeneration
 *   applied to secular equations;
 */
void
mps_context_select_algorithm (mps_context * s, mps_algorithm algorithm)
{  
  /* First set algorithm in the mps_context */
  s->algorithm = algorithm;

  switch (algorithm)
    {
    case MPS_ALGORITHM_STANDARD_MPSOLVE:
      s->mpsolve_ptr = MPS_MPSOLVE_PTR (mps_standard_mpsolve);
      break;

    case MPS_ALGORITHM_SECULAR_GA:
      s->mpsolve_ptr = MPS_MPSOLVE_PTR (mps_secular_ga_mpsolve);
      break;
    }
}


/**
 * @brief Allocate a new mps_context struct with default
 * options.
 */
mps_context *
mps_context_new ()
{
  /* Allocate the new mps_context and load default options */
  mps_context * s = (mps_context*) mps_malloc (sizeof (mps_context));
  mps_context_init (s);
  return s;
}

void
mps_context_init (mps_context * s)
{
  /* Set default streams */
  s->instr = stdin;
  s->outstr = stdout;
  s->logstr = stdout;

  /* Allocate space for the configurations */
  s->input_config  = (mps_input_configuration  *) mps_malloc (sizeof (mps_input_configuration));
  s->output_config = (mps_output_configuration *) mps_malloc (sizeof (mps_output_configuration));

  mps_set_default_values (s);

  /* Set standard precision */
  s->output_config->prec = (int) (0.9 * DBL_DIG * LOG2_10);
  MPS_DEBUG (s, "Setting prec_out to %ld digits", s->output_config->prec);

  mps_mp_set_prec (s, DBL_DIG * LOG2_10 + 1);

  s->initialized = false;
}

/**
 * @brief Free a not more useful mps_context.
 *
 * @param s the mps_context struct pointer to free.
 */
void
mps_context_free (mps_context * s)
{
  if (s->initialized)
    mps_free_data (s);

  free (s->input_config);
  free (s->output_config);

  s->active_poly = NULL;

  if (s->secular_equation)
    mps_secular_equation_free (s, MPS_POLYNOMIAL (s->secular_equation));

  /* Close input and output streams if they're not stdin, stdout and
   * stderr */
  if (s->instr != stdin && s->instr != NULL)
    fclose (s->instr);
  if (s->logstr != stderr && s->logstr != stdout && s->logstr != NULL) 
    fclose (s->logstr); 
   
   free (s);
}

void
mps_context_set_degree (mps_context * s, int n)
{
  s->deg = s->n = n;
  
  /* Check if the numer of thread is greater of the number of roots,
     and in that case decrease it */
  if (s->n_threads > s->deg)
    mps_thread_pool_set_concurrency_limit (s, s->pool, s->deg);
}

/**
 * @brief Set the monomial poly p as the input polynomial for 
 * the current equation.
 *
 * @param s The mps_context to set the monomial_poly into.
 * @param p The mps_monomial_poly to solve.
 */
void
mps_context_set_input_poly (mps_context * s, mps_polynomial * p)
{
  MPS_DEBUG_THIS_CALL;

  MPS_DEBUG (s, "Setting input poly");

  if (p->degree < 0)
  {
    mps_error (s, "Polynomial degree should be positive");
    return;
  }
  
  int i;
  s->active_poly = p;
  s->n = p->degree;

  if (!p->thread_safe)
    mps_thread_pool_set_concurrency_limit (s, s->pool, 1);

  /* Set the density or sparsity of the polynomial, if it's not
   * a user polynomial */
  if (MPS_IS_MONOMIAL_POLY (p))
    {
      mps_monomial_poly *mp = MPS_MONOMIAL_POLY (p);

      /* Deflate the polynomial if necessary */
      mps_monomial_poly_deflate (s, p);
      s->zero_roots = s->n - p->degree;
      s->n = p->degree;

      /* Check if the input polynomial is sparse or not. We can simply check if
       * the again vector is all of true values */
      p->density = MPS_DENSITY_DENSE;
      for (i = 0; i <= MPS_POLYNOMIAL (mp)->degree; ++i)
       {
         if (!mp->spar[i])
           {
             p->density = MPS_DENSITY_SPARSE;
             break;
           }
       }
    }

  mps_context_set_degree (s, p->degree);
}

/**
 * @brief Set active polynomial as a real floating point coefficient
 * polynomial of degree <code>n</code> with coefficient exactly
 * determined by components of vector coeff.
 *
 * Precisely, if \f${\mathrm coeff}\f$ is a vector of \f$n+1\f$ components,
 * \f[
 *   p(x) = \sum_{i = 0}^{n} {\mathrm coeff}_i  x^i
 * \f]
 */
int
mps_context_set_poly_d (mps_context * s, cplx_t * coeff, long unsigned int n)
{

  int i;

  /* Allocate space for a polynomial of degree n */
  mps_monomial_poly * p = mps_monomial_poly_new (s, n);

  /* Fill polynomial coefficients */
  for (i = 0; i <= n; i++)
    {
      mps_monomial_poly_set_coefficient_d (s, p, i, cplx_Re (coeff[i]),
                                           cplx_Im (coeff[i]));
    }

  mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));
  
  return 0;
}

/**
 * @brief Set active polynomial as a integer coefficient
 * polynomial of degree <code>n</code> with coefficient exactly
 * determined by components of vector coeff.
 *
 * Precisely, if \f${\mathrm coeff}\f$ is a vector of \f$n+1\f$ components,
 * \f[
 *   p(x) = \sum_{i = 0}^{n} {\mathrm coeff}_i  x^i
 * \f]
 */
int
mps_context_set_poly_i (mps_context * s, int *coeff, long unsigned int n)
{

  int i;

  /* Allocate data in mps_context to hold the polynomial of degree n */
  mps_monomial_poly * p = mps_monomial_poly_new (s, n);

  /* Fill polynomial */
  for (i = 0; i <= n; i++)
    {
      mpq_set_si (p->initial_mqp_r[i], coeff[i], 1U);
    }

  mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));

  return 0;
}

/**
 * @brief Set <code>roots[i]</code> to the i-th root of the polynomial
 * and (if it is not <code>NULL</code>) <code>radius[i]</code>
 * to the i-th inclusion radius.
 *
 * @param s The current mps_context.
 * @param roots A pointer to an array of cplx_t variables. if *roots == NULL, 
 * MPSolve will take care of allocating these for you. You are in charge to free
 * them when you don't need them anymore. 
 *
 * @param radius A pointer to an array of double where MPSolve should store the
 * inclusion radii. If *radius == NULL MPSolve will allocate those radii for you. 
 * If radius == NULL no radii will be returned. 
 */
int
mps_context_get_roots_d (mps_context * s, cplx_t ** roots, double **radius)
{
  int i;

  if (*roots == NULL)
    *roots = cplx_valloc (s->n);

  if (radius && !*radius)
    *radius = double_valloc (s->n);

  for (i = 0; i < s->n; i++)
    {
      if (radius && *radius != NULL)
        {
          if (s->lastphase == float_phase || s->lastphase == dpe_phase)
            {
              (*radius)[i] = s->root[i]->frad;
            }
          else
            {
              (*radius)[i] = rdpe_get_d (s->root[i]->drad);
            }
        }

      if (s->lastphase == mp_phase)
        {
          mpc_get_cplx ((*roots)[i], s->root[i]->mvalue);
        }
      else if (s->lastphase == float_phase)
        {
          cplx_set ((*roots)[i], s->root[i]->fvalue);
        }
      else if (s->lastphase == dpe_phase)
        {
          cdpe_get_x ((*roots)[i], s->root[i]->dvalue);
        }
    }
  return 0;
}

/**
 * @brief Get the roots computed as multiprecision complex numbers.
 *
 * @param roots A pointer to an array of mpc_t variables. if *roots == NULL, 
 * MPSolve will take care of allocate and init those for you. You are in charge to free
 * and clear them when you don't need them anymore. 
 *
 * @param radius A pointer to an array of rdpe_t where MPSolve should store the
 * inclusion radii. If *radius == NULL MPSolve will allocate those radii for you. 
 * If radius == NULL no radii will be returned. 
 */
int
mps_context_get_roots_m (mps_context * s, mpc_t ** roots, rdpe_t ** radius)
{
  int i;

  if (!*roots)
    {
      *roots = mpc_valloc (s->n);
      mpc_vinit2 (*roots, s->n, 0);
    }

  if (radius && !*radius)
    {
      *radius = rdpe_valloc (s->n);
    }

  {
    mpc_t * local_roots = *roots;
    rdpe_t * local_radius = radius ? *radius : NULL;
    
    for (i = 0; i < s->n; i++)
      {
        mpc_set_prec (local_roots[i], mpc_get_prec (s->root[i]->mvalue));
        mpc_set (local_roots[i], s->root[i]->mvalue);

        if (radius)
          rdpe_set (local_radius[i], s->root[i]->drad);
      }
  }

  return 0;
}


/**
 * @brief Set the output precision for the roots. 
 * 
 * This has differente meaning based on the output goal. 
 * If the goal is <code>MPS_OUTPUT_GOAL_ISOLATE</code>, this
 * is the maximum precision used to try to isolate the roots, 
 * but roots won't be approximated at this precision if they
 * are isolated with less precision.
 *
 * If the goal is <code>MPS_OUTPUT_GOAL_APPROXIMATE</code>, 
 * this is the minimum precision required for the roots in
 * output.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param prec The desired output precision.
 */
void 
mps_context_set_output_prec (mps_context * s, long int prec)
{
  s->output_config->prec = prec;
  rdpe_set_2dl (s->eps_out, 1.0, -prec);
}

/**
 * @brief Set the bits of precision of the input coefficients. 
 * 
 * This has meaning only for fp coefficients, and the special
 * value 0 means infinite precision.
 *
 * @param s The mps_context of the current computation.
 * @param prec The precisision to be set.
 */
void
mps_context_set_input_prec (mps_context * s, long int prec)
{
  if (!s->active_poly)
    return;
  s->active_poly->prec = prec;
}

/**
 * @brief Set the desired output format that will be used when
 * calling mps_output(). 
 *
 * @param s The <code>mps_context</code> of the current computation.
 * @param format The format chosen for the output.
 */
void 
mps_context_set_output_format (mps_context * s, mps_output_format format)
{
  s->output_config->format = format;
}

/**
 * @brief Set the output goal for the computation.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param goal The goal that will be reached before stopping.
 */
void 
mps_context_set_output_goal (mps_context * s, mps_output_goal goal)
{
  s->output_config->goal = goal;
}

/**
 * @brief Set the value of the jacobi iterations switch in the MPSolve context.
 *
 * If jacobi_iterations is true then the Ehrlich-Aberth iterations will be carried
 * out in a Jacobi fashion, otherwise Gauss-Seidel style will be employed. 
 *
 * @param s The mps_context where the value will be set
 * @param jacobi_iterations The desired value for the jacobi_iterations switch.
 */
void 
mps_context_set_jacobi_iterations (mps_context * s, mps_boolean jacobi_iterations)
{
  s->jacobi_iterations = jacobi_iterations;
}


/**
 * @brief Set the debug level in MPSolve.
 *
 * @param s The <code>mps_context</code> of the current computation.
 * @param level The debug_level to set.
 */
void 
mps_context_set_debug_level (mps_context * s, mps_debug_level level)
{
  s->debug_level = level;
  if (level)
    {
      s->DOLOG = true;
      if (!s->logstr)
        s->logstr = stderr;
    }
}

/**
 * @brief Add another debug domain to the ones displayed. This will
 * enable debug if disabled and show message from the given region.
 *
 * @param s The <code>mps_context</code> of the current computation.
 * @param level The domain to add to the already set debug_level.
 */
void 
mps_context_add_debug_domain (mps_context * s, mps_debug_level level)
{
  mps_context_set_debug_level (s, s->debug_level | level);
}

/**
 * @brief Get a pointer to the input config stored in the 
 * given mps_context.
 *
 * @param s The <code>mps_context</code> of the current computation.
 */
mps_input_configuration * 
mps_context_get_input_config (mps_context * s)
{
  return s->input_config;
}

/**
 * @brief Get a pointer to the output config stored in the
 * given mps_context.
 *
 * @param s The <code>mps_context</code> of the current computation.
 */
mps_output_configuration * 
mps_context_get_output_config (mps_context * s)
{
  return s->output_config;
}


/**
 * @brief Set logstr as the default output for logging.
 *
 * @param s The <code>mps_context</code> of the current computation.
 * @param logstr The desired stream to be used for logging. 
 */
void
mps_context_set_log_stream (mps_context * s, FILE * logstr)
{
  s->logstr = logstr;
}

/**
 * @brief Set the phase from which the computation should start.
 *
 * @param s The <code>mps_context</code> of the current computation.
 * @param phase The phase which should be chosen at the start of the
 * computation.
 */
void 
mps_context_set_starting_phase (mps_context * s, mps_phase phase)
{
  s->input_config->starting_phase = phase;
}

/**
 * @brief Get the number of zero roots in the output. Must be called
 * after mps_mpsolve() has complete.
 *
 * @param s The <code>mps_context</code> of the current computation.
 */
int 
mps_context_get_zero_roots (mps_context * s)
{
  return s->zero_roots;
}

/**
 * @brief Return true of the computation has passed the maximum
 * admitted precision, and so was unable to reach desired output
 * precision without any further information on the input
 * coefficients.
 *
 * @param s The <code>mps_context</code> of the current computation.
 */
mps_boolean 
mps_context_get_over_max (mps_context * s)
{
  return s->over_max;
}

/**
 * @brief Return true if mpsolve has encountered an error
 * in the computation.
 *
 * @param s The <code>mps_context</code> of the current computation.
 */
mps_boolean mps_context_has_errors (mps_context * s)
{
  return s->error_state;
}

/**
 * @brief Return a copy of the string describing the error msg, or
 * NULL if no error has been encountered. 
 *
 * This char array should be freed by the user.
 *
 * @param s The <code>mps_context</code> of the current computation.
 */
char * mps_context_error_msg (mps_context * s)
{
  if (s->last_error)
    return strdup (s->last_error);
  else
    return NULL;
}

/**
 * @brief Retrieve a pointer to the active polynomial being solved. 
 *
 * @return A pointer to the requested active polynomial. 
 */
mps_polynomial* 
mps_context_get_active_poly (mps_context * ctx)
{
  return ctx->active_poly;
}

/**
 * @brief Retrieve the status of the root in position i. 
 * 
 * This method can be used to obtain more insight on the status of
 * the approximations previously obtained by a call to mps_context_get_roots_m()
 * or mps_context_get_roots_d().
 *
 * @return A copy to the mps_root_status of the approximation. 
 */
mps_root_status 
mps_context_get_root_status (mps_context * ctx, int i)
{
  return ctx->root[i]->status;
}
