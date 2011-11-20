/*
 * mps_interface.c
 *
 *  Created on: 05/apr/2011
 *      Author: Leonardo Robol <robol@poisson.phc.unipi.it>
 */


#include <mps/core.h>
#include <mps/link.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/**
 * @brief Call the real polynomial (or secular equation, or whatever) solver
 * and do the computation.
 *
 * The algorithm used must be selected before this call with <code>mps_select_algorithm</code>
 * and the data (the coefficients, or whatever the algorithm may require) should be provided
 * after that.
 *
 * Roots can then be obtained with the functions <code>mps_status_get_roots_*</code>
 *
 */
void
mps_mpsolve (mps_status * s)
{
  (*s->mpsolve_ptr) (s);
}

/**
 * @brief Allocator for memory to be used in mpsolve.
 */
void *
mps_malloc (size_t size)
{
  /* fprintf (stderr, "Allocating %lu bytes of memory\n", size); */
  register void *value = malloc (size);
  if (value == 0)
    {
      fprintf (stderr, "virtual memory exhausted");
      exit (1);
    }
  return value;
}

/**
 * @brief Allocate size bytes on the stack
 */
void *
mps_alloca (size_t size)
{
  register void *value = alloca (size);
  return value;
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
mps_status_select_algorithm (mps_status * s, mps_algorithm algorithm)
{
  /* First set algorithm in the mps_status */
  s->algorithm = algorithm;

  switch (algorithm)
    {
    case MPS_ALGORITHM_STANDARD_MPSOLVE:
      s->mpsolve_ptr = MPS_MPSOLVE_PTR (mps_standard_mpsolve);
      break;

    case MPS_ALGORITHM_SECULAR_MPSOLVE:
      s->data_type = strdup ("uri");

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

      /* Check if the secular equation is allocate or if only the
       * polynomial is present. In the last case, allocate an empty
       * secular equation to hold the data during the computation. */
      if (!s->secular_equation)
	s->secular_equation = mps_secular_equation_new_raw (s, s->monomial_poly->n);

      break;
    }
}


/**
 * @brief Allocate a new mps_status struct with default
 * options.
 */
mps_status *
mps_status_new ()
{
  /* Allocate the new mps_status and load default options */
  mps_status * s = (mps_status*) mps_malloc (sizeof (mps_status));
  mps_status_init (s);
  return s;
}

void
mps_status_init (mps_status * s)
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
}

/**
 * @brief Free a not more useful mps_status.
 *
 * @param s the mps_status struct pointer to free.
 */
void
mps_status_free (mps_status * s)
{
  mps_free_data (s);

  free (s->input_config);
  free (s->output_config);

  /* Check if secular equation or monomial poly need to be freed */
  if (s->monomial_poly)
    mps_monomial_poly_free (s, s->monomial_poly);
  if (s->secular_equation)
    mps_secular_equation_free (s->secular_equation);

  free (s->data_type);

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
 * // Set a polynomial of degree n with associated mps_status* s
 * // and use the provided routines to compute newton corrections.
 * mps_status_set_poly_u(s, n,
 *   MPS_FNEWTON_PTR(mps_secular_fnewton),
 *	 MPS_DNEWTON_PTR(mps_secular_dnewton),
 *	 MPS_MNEWTON_PTR(mps_secular_mnewton));
 * @endcode
 *
 * @param s The <code>mps_status</code> struct;
 * @param n The degree of the polynomial;
 * @param fnewton The routine that performs the computation of the newton correction
 *   in floating point. It must be of the type
 *   <code>(void*)(mps_status* s, cplx_t x, double *rad, cplx_t corr, mps_boolean * again)</code>
 *   and can be passed to the function with the right casting using the macro
 *   <code>MPS_FNEWTON_PTR</code>.
 * @param dnewton The routine that performs the computation of the newton correction in
 *   <code>dpe</code> precision. It must be of the type
 *   <code>(void*)(mps_status* s, cdpe_t x, rdpe_t rad, cdpe_t corr, mps_boolean * again)</code>
 *   and can be passed to the function with the right casting using the macro
 *   <code>MPS_DNEWTON_PTR</code>.
 * @param mnewton The routine that performs the computation of the newton correction in
 *   multiprecision. It must be of the type
 *   <code>(void*)(mps_status* s, mpc_t x, rdpe_t rad, mpc_t corr, mps_boolean * again)</code>
 *   and can be passed to the function with the right casting using the macro
 *   <code>MPS_MNEWTON_PTR</code>. 
 */
int
mps_status_set_poly_u (mps_status * s, int n, mps_fnewton_ptr fnewton,
                       mps_dnewton_ptr dnewton, mps_mnewton_ptr mnewton)
{
  mps_monomial_poly *p = mps_monomial_poly_new (s, n);
  s->monomial_poly = p;

  /* Set degree and allocate data */
  mps_status_set_degree (s, n);

  /* TODO: Apart from u, what should be set here? */
  s->data_type = strdup ("uri");

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
mps_status_set_degree (mps_status * s, int n)
{
  s->deg = s->n = n;
  
  /* Check if the numer of thread is greater of the number of roots,
     and in that case decrease it */
  if (s->n_threads > s->deg)
    s->n_threads = s->deg;
}

/**
 * @brief Set the monomial poly p as the input polynomial for 
 * the current equation.
 *
 * @param s The mps_status to set the monomial_poly into.
 * @param p The mps_monomial_poly to solve.
 * @param structure The algebraic structure of the polynomial. This can
 * be, for example, <code>MPS_STRUCTURE_REAL_INTEGER</code> or similar values. 
 * What is set here will determine the fields of the poly that will be looked for data.
 */
void
mps_status_set_input_poly (mps_status * s, mps_monomial_poly * p, mps_structure structure)
{
  int i;
  s->monomial_poly = p;
  mps_status_set_degree (s, p->n);

  /* Set the right flag for the input */
  s->input_config->representation = MPS_REPRESENTATION_MONOMIAL;

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

  /* Set the mps_structure passed as input */
  s->input_config->structure = structure;
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
mps_status_set_poly_d (mps_status * s, cplx_t * coeff, long unsigned int n)
{

  int i;

  /* Allocate space for a polynomial of degree n */
  mps_monomial_poly * p = mps_monomial_poly_new (s, n);

  /* Set type to a dense, real, floating point polynomial */
  s->data_type = strdup ("dcf");

  /* Fill polynomial coefficients */
  for (i = 0; i <= n; i++)
    {
      mpc_set_cplx (p->mfpc[i], coeff[i]);
      cplx_set (p->fpc[i], coeff[i]);
      cdpe_set_x (p->dpc[i], coeff[i]);
    }
  
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
mps_status_set_poly_i (mps_status * s, int *coeff, long unsigned int n)
{

  int i;

  /* Allocate data in mps_status to hold the polynomial of degree n */
  mps_monomial_poly * p = mps_monomial_poly_new (s, n);

  /* Dense, real, integer coefficients */
  s->data_type = strdup ("drq");

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
mps_status_get_roots_d (mps_status * s, cplx_t * roots, double *radius)
{
  int i;
  for (i = 0; i < s->n; i++)
    {

      if (radius != NULL)
        {
          if (s->lastphase == float_phase || s->lastphase == dpe_phase)
            {
              radius[i] = s->frad[i];
            }
          else
            {
              radius[i] = rdpe_get_d (s->drad[i]);
            }

        }

      if (s->lastphase == mp_phase)
        {
          mpc_get_cplx (roots[i], s->mroot[i]);
        }
      else if (s->lastphase == float_phase)
        {
          cplx_set (roots[i], s->froot[i]);
        }
      else if (s->lastphase == dpe_phase)
        {
          cdpe_get_x (roots[i], s->droot[i]);
        }
    }
  return 0;
}

/**
 * @brief Get the roots computed as multiprecision complex numbers.
 */
int
mps_status_get_roots_m (mps_status * s, mpc_t * roots, rdpe_t * radius)
{
  /* TODO: Implement mps_get_roots_m() */
  return 0;
}
