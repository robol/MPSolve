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
  /* Select ga if we are in the case of monomial input, since is the only algorithm
   * that we can use. */
  if (MPS_INPUT_CONFIG_IS_MONOMIAL (mps_context_get_input_config (s)) && algorithm == MPS_ALGORITHM_SECULAR_MPSOLVE)
    {
      algorithm = MPS_ALGORITHM_SECULAR_GA;
      MPS_DEBUG_WITH_INFO (s, "Selecting algorithm MPS_ALGORITHM_SECULAR_GA since MPS_ALGORITHM_SECULAR_MPSOLVE is not available for monomial input");
    }
  
  /* First set algorithm in the mps_context */
  s->algorithm = algorithm;

  switch (algorithm)
    {
    case MPS_ALGORITHM_STANDARD_MPSOLVE:
      s->mpsolve_ptr = MPS_MPSOLVE_PTR (mps_standard_mpsolve);
      break;

    case MPS_ALGORITHM_SECULAR_MPSOLVE:

      /* Set custom routines for newton quotient computation */
      s->fnewton_usr = MPS_FNEWTON_PTR (mps_secular_fnewton);
      s->dnewton_usr = MPS_DNEWTON_PTR (mps_secular_dnewton);
      s->mnewton_usr = MPS_MNEWTON_PTR (mps_secular_mnewton);

      /* Set other custom functions */
      s->check_data_usr = MPS_CHECK_DATA_PTR (mps_secular_check_data);

      /* Set starting point custom routine */
      s->fstart_usr = MPS_FSTART_PTR (mps_secular_fstart);
      s->dstart_usr = MPS_DSTART_PTR (mps_secular_dstart);

      /* Set default routine for mpsolve loop */
      s->mpsolve_ptr = MPS_MPSOLVE_PTR (mps_standard_mpsolve);

      break;

    case MPS_ALGORITHM_SECULAR_GA:
      
      /* Nothing to be done for now, other than selecting the correct
       * loop. */
      s->mpsolve_ptr = MPS_MPSOLVE_PTR (mps_secular_ga_mpsolve);

      /* Set custom routines for newton quotient computation, that will
       * become useful in mps_improve () */
      s->fnewton_usr = MPS_FNEWTON_PTR (mps_secular_fnewton);
      s->dnewton_usr = MPS_DNEWTON_PTR (mps_secular_dnewton);
      s->mnewton_usr = MPS_MNEWTON_PTR (mps_secular_mnewton);

      rdpe_set_2dl (s->eps_out, 1.0, -s->output_config->prec);

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
  s->input_config->prec = 0;

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

  /* Check if secular equation or monomial poly need to be freed */
  if (s->monomial_poly)
    mps_monomial_poly_free (s, s->monomial_poly);
  if (s->secular_equation)
    mps_secular_equation_free (s->secular_equation);

  /* Close input and output streams if they're not stdin, stdout and
   * stderr */
  if (s->instr != stdin && s->instr != NULL)
    fclose (s->instr);
  if (s->logstr != stderr && s->logstr != stdout && s->logstr != NULL) 
    fclose (s->logstr); 
   
   free (s);
}

/**
 * @brief Set active poly as a user poly, providing routines to compute
 * newton corrections.
 *
 * This is an example of call to this function:
 * @code
 * // Set a polynomial of degree n with associated mps_context* s
 * // and use the provided routines to compute newton corrections.
 * mps_context_set_poly_u(s, n,
 *   MPS_FNEWTON_PTR(mps_secular_fnewton),
 *	 MPS_DNEWTON_PTR(mps_secular_dnewton),
 *	 MPS_MNEWTON_PTR(mps_secular_mnewton));
 * @endcode
 *
 * @param s The <code>mps_context</code> struct;
 * @param n The degree of the polynomial;
 * @param fnewton The routine that performs the computation of the newton correction
 *   in floating point. It must be of the type
 *   <code>(void*)(mps_context* s, cplx_t x, double *rad, cplx_t corr, mps_boolean * again)</code>
 *   and can be passed to the function with the right casting using the macro
 *   <code>MPS_FNEWTON_PTR</code>.
 * @param dnewton The routine that performs the computation of the newton correction in
 *   <code>dpe</code> precision. It must be of the type
 *   <code>(void*)(mps_context* s, cdpe_t x, rdpe_t rad, cdpe_t corr, mps_boolean * again)</code>
 *   and can be passed to the function with the right casting using the macro
 *   <code>MPS_DNEWTON_PTR</code>.
 * @param mnewton The routine that performs the computation of the newton correction in
 *   multiprecision. It must be of the type
 *   <code>(void*)(mps_context* s, mpc_t x, rdpe_t rad, mpc_t corr, mps_boolean * again)</code>
 *   and can be passed to the function with the right casting using the macro
 *   <code>MPS_MNEWTON_PTR</code>. 
 */
int
mps_context_set_poly_u (mps_context * s, int n, mps_fnewton_ptr fnewton,
                       mps_dnewton_ptr dnewton, mps_mnewton_ptr mnewton)
{
  mps_monomial_poly *p = mps_monomial_poly_new (s, n);
  s->monomial_poly = p;

  /* Set degree and allocate data */
  mps_context_set_degree (s, n);

  /* Set functions */
  s->fnewton_usr = fnewton;
  s->dnewton_usr = dnewton;
  s->mnewton_usr = mnewton;

  s->input_config->structure = MPS_STRUCTURE_REAL_INTEGER;
  s->input_config->density = MPS_DENSITY_USER;
  s->input_config->representation = MPS_REPRESENTATION_MONOMIAL;

  return 0;
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
mps_context_set_input_poly (mps_context * s, mps_monomial_poly * p)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  s->monomial_poly = p;
  mps_context_set_degree (s, p->n);

  /* Set the right flag for the input */
  s->input_config->representation = MPS_REPRESENTATION_MONOMIAL;

  /* Set the mps_structure passed as input */
  s->input_config->structure = p->structure;

  /* Set the density or sparsity of the polynomial, if it's not
   * a user polynomial */
  if (!MPS_INPUT_CONFIG_IS_USER (s->input_config))
    {
      /* Check if the input polynomial is sparse or not. We can simply check if
       * the again vector is all of true values */
      s->input_config->density = MPS_DENSITY_DENSE;
      for (i = 0; i <= p->n; ++i)
	{
	  if (!p->spar[i])
	    {
	      s->input_config->density = MPS_DENSITY_SPARSE;
	      break;
	    }
	}
    }
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

  mps_context_set_input_poly (s, p);
  
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

  return 0;
}

/**
 * @brief Set <code>roots[i]</code> to the i-th root of the polynomial
 * and (if it is not <code>NULL</code>) <code>radius[i]</code>
 * to the i-th inclusion radius.
 */
int
mps_context_get_roots_d (mps_context * s, cplx_t * roots, double *radius)
{
  int i;
  for (i = 0; i < s->n; i++)
    {

      if (radius != NULL)
        {
          if (s->lastphase == float_phase || s->lastphase == dpe_phase)
            {
              radius[i] = s->root[i]->frad;
            }
          else
            {
              radius[i] = rdpe_get_d (s->root[i]->drad);
            }

        }

      if (s->lastphase == mp_phase)
        {
          mpc_get_cplx (roots[i], s->root[i]->mvalue);
        }
      else if (s->lastphase == float_phase)
        {
          cplx_set (roots[i], s->root[i]->fvalue);
        }
      else if (s->lastphase == dpe_phase)
        {
          cdpe_get_x (roots[i], s->root[i]->dvalue);
        }
    }
  return 0;
}

/**
 * @brief Get the roots computed as multiprecision complex numbers.
 */
int
mps_context_get_roots_m (mps_context * s, mpc_t * roots, rdpe_t * radius)
{
  int i;

  mps_copy_roots (s);

  for (i = 0; i < s->n; i++)
    {
      mpc_set_prec (roots[i], mpc_get_prec (s->root[i]->mvalue));
      mpc_set (roots[i], s->root[i]->mvalue);
      rdpe_set (radius[i], s->root[i]->drad);
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
  s->input_config->prec = prec;
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
