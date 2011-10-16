/*
 * mps_interface.c
 *
 *  Created on: 05/apr/2011
 *      Author: Leonardo Robol <robol@poisson.phc.unipi.it>
 */

/**
 * 	@file
 * 	@brief Implementation of the routines to interact with MPSolve
 * 	as a library.
 */

/**
 * @mainpage Using MPSolve as a library
 *
 * @section Installation Installing MPSolve system-wide
 * First, you need to get MPSolve. You can get the latest release via <code>git</code>
 * or download it via <code>http</code> grabbing it at http://www.dm.unipi.it/...
 * If you downloaded the source tarball this operation is pretty straightforward.
 * You can simply unpack it and then
 * @code
 *   make
 *   [sudo] make install
 * @endcode
 *
 * These commands will install the library <code>libmps.so</code> in your system library
 * directory. In this way you will be able compile your source file using a command similar
 * to
 * @code
 * gcc -o myprogram -lmps -lgmp -lm myprogram.c.
 * @endcode
 *
 * @section Interface Using the libmps interface
 *
 * The library provides some useful routine to interact with the polynomial solver. Most of
 * them are designed to handle polynomial definition and are implemented in mps_interface.c
 *
 */

#include <mps/core.h>
#include <mps/poly.h>
#include <mps/link.h>
#include <mps/secular.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>

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
  int i;

  /* First set algorithm in the mps_status */
  s->algorithm = algorithm;

  switch (algorithm)
    {
    case MPS_ALGORITHM_STANDARD_MPSOLVE:
      s->mpsolve_ptr = MPS_MPSOLVE_PTR (mps_standard_mpsolve);
      break;

    case MPS_ALGORITHM_SECULAR_MPSOLVE:
      s->data_type = "uri";

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

      rdpe_set_2dl (s->eps_out, 1.0, -s->prec_out * LOG2_10);

      break;
    }
}

/**
 * @brief Allocate polynomial related variables directly in mps_status.
 */
void
mps_status_allocate_poly_inplace (mps_status * s, int n)
{

  int i;

  /* If n is provided, then we should allocate variables for a polynomial
   * of degree n. If it is <= 0 then we can assume that is already set
   * to the right vaule.
   */
  if (n > 0)
    {
      s->deg = s->n = n;
    }

  MPS_DEBUG (s, "Allocating polynomial in place");

  if (!s->data_type)
    s->data_type = (char *) malloc (sizeof (char) * 3);

  s->spar = mps_boolean_valloc (s->deg + 2);

  s->fpr = double_valloc (s->deg + 1);
  s->fpc = cplx_valloc (s->deg + 1);

  s->dpr = rdpe_valloc (s->deg + 1);
  s->dpc = cdpe_valloc (s->deg + 1);

  s->mip_r = mpz_valloc (s->deg + 1);
  s->mip_i = mpz_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    {
      mpz_init (s->mip_r[i]);
      mpz_init (s->mip_i[i]);
    }

  s->mqp_r = mpq_valloc (s->deg + 1);
  s->mqp_i = mpq_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    {
      mpq_init (s->mqp_r[i]);
      mpq_init (s->mqp_i[i]);
    }

  s->mfpr = mpf_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    mpf_init2 (s->mfpr[i], s->prec_in);
  s->mfpc = mpc_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    mpc_init2 (s->mfpc[i], s->prec_in);

  /* Create the status and set all the roots as uncertain */
  s->status = (char (*)[3]) char_valloc (3 * s->deg);
  for (i = 0; i < s->deg; i++)
    {
      s->status[i][1] = 'w';
      s->status[i][2] = 'u';
    }

}

void
mps_status_free_poly_inplace (mps_status * s)
{
  free (s->spar);
  cplx_vfree (s->fpc);
  rdpe_vfree (s->dpr);
  rdpe_vfree (s->dpc);

  mpz_vclear (s->mip_r, s->n);
  mpz_vclear (s->mip_i, s->n);

  mpz_vfree (s->mip_r);
  mpz_vfree (s->mip_i);

  mpq_vclear (s->mqp_r, s->n);
  mpq_vclear (s->mqp_i, s->n);

  mpq_vfree (s->mqp_r);
  mpq_vfree (s->mqp_i);

  mpf_vclear (s->mfpr, s->n);
  mpc_vclear (s->mfpc, s->n);

  mpf_vfree (s->mfpr);
  mpf_vfree (s->mfpc);
}

/**
 * @brief Allocate a new mps_status struct with default
 * options.
 */
mps_status *
mps_status_new ()
{
  /* Allocate the new mps_status and load default options */
  mps_status *s = (mps_status *) malloc (sizeof (mps_status));
  mps_set_default_values (s);

  /* Set default streams */
  s->instr = stdin;
  s->outstr = stdout;
  s->logstr = stdout;

  /* Set standard precision */
  s->prec_out = (int) (0.9 * DBL_DIG * LOG2_10);
  MPS_DEBUG (s, "Setting prec_out to %d digits", s->prec_out);
  s->prec_in = 0;

  return s;
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

  /* Set degree and allocate data */
  mps_status_set_degree (s, n);
  mps_allocate_data (s);

  /* TODO: Apart from u, what should be set here? */
  s->data_type = "uri";

  /* Set functions */
  s->fnewton_usr = fnewton;
  s->dnewton_usr = dnewton;
  s->mnewton_usr = mnewton;

  return 0;
}

void
mps_status_set_degree (mps_status * s, int n)
{
  s->deg = s->n = n;
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
  mps_status_allocate_poly_inplace (s, n);

  /* Set type to a dense, real, floating point polynomial */
  s->data_type[0] = 'd';
  s->data_type[1] = 'c';
  s->data_type[2] = 'f';

  /* Fill polynomial coefficients */
  for (i = 0; i <= n; i++)
    {
      mpc_set_cplx (s->mfpc[i], coeff[i]);
    }

  /* Allocate space for computation related data */
  mps_allocate_data (s);

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
  mps_status_allocate_poly_inplace (s, n);

  /* Dense, real, integer coefficients */
  s->data_type[0] = 'd';
  s->data_type[1] = 'r';
  s->data_type[2] = 'i';

  /* Fill polynomial */
  for (i = 0; i <= n; i++)
    {
      mpz_set_si (s->mip_r[i], coeff[i]);
    }

  /* Allocate data for the computation */
  mps_allocate_data (s);

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
  /* TODO: Implement mps_get_roots_d() */
  return 0;
}
