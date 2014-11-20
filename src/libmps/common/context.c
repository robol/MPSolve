/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#define _GNU_SOURCE

#include <mps/mps.h>
#include <string.h>

static mps_context ** context_factory = NULL;
static int context_factory_size = 0;
static pthread_mutex_t context_factory_mutex = PTHREAD_MUTEX_INITIALIZER;

#define MPS_CONTEXT_FACTORY_MAXIMUM_SIZE 0

long int
mps_context_get_minimum_precision (mps_context * s)
{
  return s->minimum_gmp_precision;
}

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
 * @brief Select the starting strategy used to dispose the initial approximations. 
 */
void
mps_context_select_starting_strategy (mps_context * s, mps_starting_strategy strategy)
{
  s->starting_strategy = strategy;
}

static void
mps_context_init (mps_context * s)
{
  mpf_t test;

  /* Set default streams */
  s->instr = stdin;
  s->outstr = stdout;
  s->logstr = stdout;

  /* Allocate space for the configurations */
  s->input_config = (mps_input_configuration*)mps_malloc (sizeof(mps_input_configuration));
  s->output_config = (mps_output_configuration*)mps_malloc (sizeof(mps_output_configuration));

  mps_set_default_values (s);

  /* Find minimum GMP supported precision */
  mpf_init2 (test, 1);
  s->minimum_gmp_precision = mpf_get_prec (test);
  mpf_clear (test);

  /* Set standard precision */
  s->output_config->prec = (int)(0.9 * DBL_DIG * LOG2_10);
  MPS_DEBUG (s, "Setting prec_out to %ld digits", s->output_config->prec);

  mps_mp_set_prec (s, DBL_DIG * LOG2_10 + 1);

  s->initialized = false;
  s->exit_required = false;
}

/**
 * @brief Allocate a new mps_context struct with default
 * options.
 */
mps_context *
mps_context_new ()
{
  mps_context * ctx = NULL;

  /* Try to get a low cost context from the factory */
  pthread_mutex_lock (&context_factory_mutex);
  if (context_factory_size > 0)
    {
      /* Pop out a context */
      ctx = context_factory[--context_factory_size];

      if (context_factory_size)
        context_factory = mps_realloc (context_factory,
                                       sizeof(mps_context*) * context_factory_size);
      else
        {
          free (context_factory);
          context_factory = NULL;
        }
    }
  pthread_mutex_unlock (&context_factory_mutex);

  /* Allocate the new mps_context and load default options */
  if (!ctx)
    {
      ctx = (mps_context*)mps_malloc (sizeof(mps_context));
      mps_context_init (ctx);
    }

  return ctx;
}


/**
 * @brief Free a not more useful mps_context.
 *
 * @param s the mps_context struct pointer to free.
 */
void
mps_context_free (mps_context * s)
{
  /* Close input and output streams if they're not stdin, stdout and
   * stderr. For the case in which this context will re-used, set them
   * to their default values. */
  if (s->instr != stdin && s->instr != NULL)
    fclose (s->instr);
  if (s->logstr != stderr && s->logstr != stdout && s->logstr != NULL)
    fclose (s->logstr);

  s->instr = stdin;
  s->logstr = stderr;

  /* There's no need to resize bmpc since they will be allocated on demand.
   * We free them here to correct bad assumptions on the size of this
   * vector. */
  free (s->bmpc);
  s->bmpc = NULL;

  pthread_mutex_lock (&context_factory_mutex);

  if (context_factory_size < MPS_CONTEXT_FACTORY_MAXIMUM_SIZE)
    {
      context_factory = mps_realloc (context_factory,
                                     sizeof(mps_context*) * (context_factory_size + 1));
      context_factory[context_factory_size++] = s;
      pthread_mutex_unlock (&context_factory_mutex);
      return;
    }
  pthread_mutex_unlock (&context_factory_mutex);

  if (s->initialized)
    mps_free_data (s);

  mps_thread_pool_free (s, s->pool);

  free (s->input_config);
  free (s->output_config);

  s->active_poly = NULL;

  if (s->secular_equation)
    mps_secular_equation_free (s, MPS_POLYNOMIAL (s->secular_equation));

  free (s);
}

void
mps_context_abort (mps_context * s)
{
  s->exit_required = true;
}

static void
mps_context_shrink (mps_context * s, int n)
{
  int i;

  for (i = n; i < s->n - s->zero_roots; i++)
    {
      mps_approximation_free (s, s->root[i]);
    }

  s->root = mps_realloc (s->root, sizeof(mps_approximation*) * n);

  s->order = mps_realloc (s->order, sizeof(int) * n);

  s->fppc1 = mps_realloc (s->fppc1, sizeof(cplx_t) * (n + 1));

  for (i = n + 1; i <= s->n - s->zero_roots; i++)
    mpc_clear (s->mfpc1[i]);

  s->mfpc1 = mps_realloc (s->mfpc1, sizeof(mpc_t) * (n + 1));

  for (i = n + 1; i <= s->n- s->zero_roots; i++)
    mpc_clear (s->mfppc1[i]);

  s->mfppc1 = mps_realloc (s->mfppc1, sizeof(mpc_t) * (n + 1));

  /* temporary vectors */
  s->spar1 = mps_realloc (s->spar1, sizeof(mps_boolean) * (n + 2));
  s->again_old = mps_realloc (s->again_old, sizeof(mps_boolean) * (n));

  s->fap1 = mps_realloc (s->fap1, sizeof(double) * (n + 1));
  s->fap2 = mps_realloc (s->fap2, sizeof(double) * (n + 1));

  s->dap1 = mps_realloc (s->dap1, sizeof(rdpe_t) * (n + 1));
  s->dpc1 = mps_realloc (s->dpc1, sizeof(cdpe_t) * (n + 1));
  s->dpc2 = mps_realloc (s->dpc2, sizeof(cdpe_t) * (n + 1));

  /* Setting some default here, that were not settable because we didn't know
   * the degree of the polynomial */
  for (i = 0; i < n; i++)
    s->root[i]->wp = DBL_DIG * LOG2_10;
}

static void
mps_context_expand (mps_context * s, int n)
{
  int i;
  long int previous_prec = mpc_get_prec (s->mfpc1[0]);

  s->root = mps_realloc (s->root, sizeof(mps_approximation*) * n);
  for (i = s->n - s->zero_roots; i < n; i++)
    {
      s->root[i] = mps_approximation_new (s);
    }

  s->order = mps_realloc (s->order, sizeof(int) * n);

  s->fppc1 = mps_realloc (s->fppc1, sizeof(cplx_t) * (n + 1));
  s->mfpc1 = mps_realloc (s->mfpc1, sizeof(mpc_t) * (n + 1));

  for (i = s->n + 1 - s->zero_roots; i < n + 1; i++)
    mpc_init2 (s->mfpc1[i], previous_prec);

  s->mfppc1 = mps_realloc (s->mfppc1, sizeof(mpc_t) * (n + 1));
  for (i = s->n + 1- s->zero_roots; i <= n; i++)
    mpc_init2 (s->mfppc1[i], previous_prec);

  /* temporary vectors */
  s->spar1 = mps_realloc (s->spar1, sizeof(mps_boolean) * (n + 2));
  s->again_old = mps_realloc (s->again_old, sizeof(mps_boolean) * (n));

  s->fap1 = mps_realloc (s->fap1, sizeof(double) * (n + 1));
  s->fap2 = mps_realloc (s->fap2, sizeof(double) * (n + 1));

  s->dap1 = mps_realloc (s->dap1, sizeof(rdpe_t) * (n + 1));
  s->dpc1 = mps_realloc (s->dpc1, sizeof(cdpe_t) * (n + 1));
  s->dpc2 = mps_realloc (s->dpc2, sizeof(cdpe_t) * (n + 1));

  /* Setting some default here, that were not settable because we didn't know
   * the degree of the polynomial */
  for (i = 0; i < n; i++)
    s->root[i]->wp = DBL_DIG * LOG2_10;
}

void
mps_context_resize (mps_context * s, int n)
{
  /* We're excluding the case n == s->n that, clearly, doesn't
   * need any operation at all. */
  if (n > s->n)
    mps_context_expand (s, n);
  if (n < s->n)
    mps_context_shrink (s, n);
}

void
mps_context_set_degree (mps_context * s, int n)
{
  if (s->initialized)
    {
      if (s->secular_equation != NULL)
	{
	  mps_secular_equation_free (s, MPS_POLYNOMIAL (s->secular_equation));
	  s->secular_equation = NULL;
	}

      mps_context_resize (s, n);
    }

  s->deg = s->n = n;

  /* Check if the numer of thread is greater of the number of roots,
     and in that case decrease it */
  if (s->n_threads > s->deg)
    {
      MPS_DEBUG_WITH_INFO (s, "Adjusting concurrency limit to %d", s->deg);
      mps_thread_pool_set_concurrency_limit (s, s->pool, s->deg);
    }
  
  /* If a secular equation is present in the old context we should free it
   * now so it will be reallocated on the first call to the algorithm. */
  if (s->secular_equation && MPS_POLYNOMIAL (s->secular_equation)->degree < n)
    mps_secular_equation_free (s, MPS_POLYNOMIAL (s->secular_equation));
  s->secular_equation = NULL;

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
  MPS_DEBUG_THIS_CALL (s);

  MPS_DEBUG (s, "Setting input poly");

  if (p->degree < 0)
    {
      mps_error (s, "Polynomial degree should be positive");
      return;
    }

  int i;
  s->active_poly = p;

  if (!p->thread_safe)
    mps_thread_pool_set_concurrency_limit (s, s->pool, 1);

  /* Set the density or sparsity of the polynomial, if it's not
   * a user polynomial */
  if (MPS_IS_MONOMIAL_POLY (p))
    {
      int original_degree = p->degree;
      mps_monomial_poly *mp = MPS_MONOMIAL_POLY (p);

      /* Deflate the polynomial if necessary */
      mps_monomial_poly_deflate (s, p);
      s->zero_roots = original_degree - p->degree;

      MPS_DEBUG_WITH_INFO (s, "Degree = %d", p->degree);

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
 * @brief Retrieve a pointer to the current approximations in the context.
 *
 * @param ctx The current mps_context.
 * @return A vector of n mps_approximation pointer,
 * where n is the degree of the current polynomial
 * that can be retrieve with mps_context_get_degree().
 *
 * Note that the value returned is a copy of the original approximations, so
 * it should be freed when not needed anymore.
 *
 * Also note that mps_approximation_free() need the context where the approximations
 * were created to proceed, so you should free those approximaitions _before_
 * freeing the context.
 */
mps_approximation **
mps_context_get_approximations (mps_context * ctx)
{
  mps_approximation ** approximations = NULL;
  int i;

  if (!ctx->root)
    {
      /* This means that an error occurred in the previous computation and the
       * polynomial has not been solved (or mps_mpsolve has not been called at all). */
      return NULL;
    }
  else
    approximations = mps_newv (mps_approximation *, ctx->n + ctx->zero_roots);

  for (i = 0; i < ctx->n; i++)
    {
      approximations[i] = mps_approximation_copy (ctx, ctx->root[i]);

      /* Copy relevant data from the multiprecision values */
      mpc_get_cdpe (approximations[i]->dvalue, approximations[i]->mvalue);
      mpc_get_cplx (approximations[i]->fvalue, approximations[i]->mvalue);
      approximations[i]->frad = rdpe_get_d (approximations[i]->drad);
    }

  for (i = ctx->n; i < ctx->n + ctx->zero_roots; i++)
    {
      approximations[i] = mps_approximation_new (ctx);
      approximations[i]->status = MPS_ROOT_STATUS_APPROXIMATED;

      mpc_set_ui (approximations[i]->mvalue, 0U, 0U);
      cdpe_set (approximations[i]->dvalue, cdpe_zero);
      cplx_set (approximations[i]->fvalue, cplx_zero);

      rdpe_set (approximations[i]->drad, rdpe_zero);
      approximations[i]->frad = 0.0;
    }

  return approximations;
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
  else
    mps_polynomial_set_input_prec (s, s->active_poly, prec);
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

  /* Special handling of GNUPLOT format */
  if (format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      s->gnuplot_format = "xyerrorbars";
    }
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

/**
 * @brief Set the internal flag "avoid_multiprecision" to the specified value. 
 *
 * If avoid_multiprecision is true MPSolve will not enter a multiprecision stage
 * thus making impossible the computation of more digits than the one that are
 * representable in standard floating point. 
 *
 * This may be a useful flag to approximate roots of very ill-conditioned polynomials
 * of high degree when no strict isolation is required. In this case it's possible to 
 * obtain very good approximations that are not Newton-isolated but are still satisfactory. 
 */
void
mps_context_set_avoid_multiprecision (mps_context * s, mps_boolean avoid_multiprecision)
{
  s->avoid_multiprecision = avoid_multiprecision;
}

/**
 * @brief Enable the "crude" only approximation mode of MPSolve. 
 *
 * If this mode is activated MPSolve will only perform a basic Aberth iteration
 * in floating point and then exit. Note that the output result will still be
 * guaranteed but in general it will not be possible to reach arbitrary precision
 * and the results may be quite far from the roots for bad conditioned polynomials. 
 */
void
mps_context_set_crude_approximation_mode (mps_context * s, mps_boolean crude_approximation_mode)
{
  s->crude_approximation_mode = crude_approximation_mode;
}
