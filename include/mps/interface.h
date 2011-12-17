/* 
 * File:   interface.h
 * Author: leonardo
 *
 * Created on 23 aprile 2011, 10.35
 */

#ifndef MPS_INTERFACE_H
#define	MPS_INTERFACE_H

#ifdef	__cplusplus
extern "C"
{
#endif

  /* Octave module workardound */
#ifdef __UNDEF_CPLUSPLUS
#undef __cplusplus
#endif

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <gmp.h>
#include <mps/mps.h>

  /**
   * @file
   * @brief Simple routines used to interact with MPSolve without going into the internals.
   */

  /* Since this is the interface header it contains also the documentation that doxygen
   * will output in HTML form. */

  /**
   * @mainpage General documentation
   *
   * @section About What is MPSolve
   * MPSolve is a C library that allow to find solution to univariate polynomial equations 
   * with arbitrary precision. 
   *
   * More precisely, MPSolve can handle polynomials but not necessarly in their monomial
   * form. Support is now given even for secular equation and more implicit representation
   * are scheduled to be added later.
   *
   * @section Installation Installing MPSolve system-wide
   * First, you need to get MPSolve. You can get the latest release via <code>git</code>
   * or download it via <code>http</code> grabbing it at http://www.dm.unipi.it/...
   * If you downloaded the source tarball this operation is pretty straightforward.
   * You can simply unpack it and then
   * @code
   *   ./configure
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
   * @section ExtendedTypes MPSolve extended arithmetic types
   * To perform its computations MPSolve uses some more types than standard C does.
   * You will need to interact with these to give and obtain data from MPSolve. 
   * 
   * There are three main category of types that you will be required to deal with:
   * -# Standard floating point doubles;
   * -# DPE types;
   * -# Multiprecision GMP types
   *
   * @subsection FloatingPoint Floating point simple types
   * There is clearly nothing to explain about floating point doubles, but MPSolve need to deal
   * with complex floating point and it uses a type called <code>cplx_t</code> that is nothing
   * more than a struct with two double field, <code>r</code> and <code>i</code> that represents
   * the real and imaginary part of the given complex number. 
   * You should access these fields with macro provided in this way:
   * @code
   * cplx_t my_complex_number;
   * double theta = 0.5;
   * cplx_Re (my_complex_number) = cos (theta);
   * cplx_Im (my_complex_number) = sin (theta);
   * @endcode
   *
   * Some standard complex number are provided for convenicence, such as <code>cplx_one</code>
   * and <code>cplx_zero</code>.
   * 
   * To get an idea of all the routines that you can use you can read the <code>mt.h</code> header
   * file. But for a simple start, the routines that you will need are cplx_add(),
   * cplx_sub(), cplx_mul() and cplx_div(). I think that it's pretty clear what they do.
   * 
   * @subsection DPE The DPE types
   * The <code>DPE</code> are have the same precision of the <code>double</code>s but allow
   * the exponent to be much larger; the exponent is stored as a <code>long</code> value so
   * the maximum reachable one is <code>LONG_MAX</code>, while the minimum is <code>LONG_MIN</code>.
   *
   * In the <code>mt</code> library that is embedded in MPSolve two version of the <code>DPE</code>
   * types are provided: the real and the complex one. They are called <code>RDPE</code> and <code>CDPE</code>,
   * respectively.
   * 
   * The function used to handle these types are almost the same of the complex one case, so we will
   * not cover they extensively here.
   *
   * @subsection GMP GMP types
   * GMP types are used to represent multiprecision data that is used to get arbitrary high precision
   * approximation of the roots. See http://gmplib.org for details on the use of GMP. 
   *
   * @section Interface Using the libmps interface
   *
   * The library provides some useful routine to interact with the polynomial solver. Most of
   * them are designed to handle polynomial definition and are implemented in mps_interface.c
   *
   */

  /*
   *    ====== ROUTINES EXPOSED TO THE INTERFACE ======
   */

  /* functions in mps_defaults.c */
  void mps_set_default_values (mps_status * s);

  /* Functions in mps_main.c */
  void mps_mpsolve (mps_status * s);
  void mps_standard_mpsolve (mps_status * s);

  /* functions in mps_interface.c */
  void * mps_malloc (size_t size);
  void * mps_alloca (size_t size);

  /* Macros to init pointer and/or vectors in a convenient way */
#define mps_new(type) ((type *) mps_malloc (sizeof (type)))
#define mps_newv(type, n) ((type *) mps_malloc (sizeof (type) * n))

  mps_status * mps_status_new (void);
  void mps_status_init (mps_status * s);
  void mps_status_free (mps_status * s);
  int mps_status_set_poly_d (mps_status * s, cplx_t * coeff,
			     long unsigned int n);
  void mps_status_set_input_poly (mps_status * s, mps_monomial_poly * p, mps_structure structure);
  int mps_status_set_poly_i (mps_status * s, int *coeff, long unsigned int n);
  int mps_status_get_roots_d (mps_status * s, cplx_t * roots, double *radius);
  int mps_status_set_poly_u (mps_status * s, int n, mps_fnewton_ptr fnewton,
			     mps_dnewton_ptr dnewton,
			     mps_mnewton_ptr mnewton);
  void mps_status_allocate_poly_inplace (mps_status * s, int n);
  void mps_status_select_algorithm (mps_status * s, mps_algorithm algorithm);
  void mps_status_set_degree (mps_status * s, int n);
  int mps_status_get_roots_d (mps_status * s, cplx_t * roots, double *radius);
  int mps_status_get_roots_m (mps_status * s, mpc_t * roots, rdpe_t * radius);

#ifdef	__cplusplus
}
#endif

#ifdef __UNDEF_CPLUSPLUS
}
#endif

#endif                          /* MPS_INTERFACE_H */
