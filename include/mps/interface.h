/* 
 * File:   interface.h
 * Author: leonardo
 *
 * Created on 23 aprile 2011, 10.35
 */

#ifndef MPS_INTERFACE_H
#define MPS_INTERFACE_H

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <gmp.h>
#include <mps/mps.h>

#ifdef  __cplusplus
extern "C"
{
#endif

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
   * There is clearly nothing to explain about floating point doubles, but MPSolve needs to deal
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
   * them are designed to handle polynomial definition and are implemented in <code>interface.c</code>
   *
   * This is a typical example of how you'll be using MPSolve.
   *
   * @code
   * // Select the degree
   * int n = 4;
   *
   * // Allocate a new mps_context that hold the status of a computation.
   * // Its field should never be accessed directly but only via appropriate
   * //functions.
   * mps_context * status = mps_context_new ();
   * 
   * // Create a polynomial that will be solved
   * mps_monomial_poly * poly = mps_monomial_poly_new (status, n);
   *
   * // Set the coefficients. We will solve x^n - 1 in here
   * mps_monomial_poly_set_coefficient_int (status, poly, 0, -1, 0);
   * mps_monomial_poly_set_coefficient_int (status, poly, n,  1, 0);
   *
   * // Select some common output options, i.e. 512 bits of precision
   * // (more or less 200 digits guaranteed) and approximation goal.
   * mps_context_set_output_precision (status, 512);
   * mps_context_set_output_goal (status, MPS_OUTPUT_GOAL_APPROXIMATE);
   *
   * // Solve the polynomial
   * mps_context_set_input_poly (status, poly);
   * mps_mpsolve (status);
   * 
   * // Get the roots in a <code>cplx_t</code> vector. Please note that
   * // this make completely useless to have asked 512 bits of output
   * // precision, and you should use mps_context_get_roots_m() to get
   * // multiprecision approximation of the roots.
   * cplx_t * results = cplx_valloc (n);
   * mps_context_get_roots_d (status, &results, NULL);
   *
   * // Free the data used. This will free the monomial_poly if you have
   * // not done it by yourself.
   * mps_context_free (status);
   * cplx_vfree (results);
   * @endcode
   *
   * As pointed out in the comments, this piece of code is not that smart, since it asks
   * for a lot of digits in output, but then discard all this good by copying the roots
   * in floating point <code>cplx_t</code> types. 
   *
   * Please see the documentation of mps_context_get_roots_m() for a better way to get
   * the results. It was not used in example in order to keep it short and clear. 
   *
   * @subsection inclusion How to include libmps
   *
   * In general libmps will be usable by including the header file <code>mps/mps.h</code>.
   * Please note that almost all structures defined in MPSolve will be available only as
   * incomplete declarations, so they should be used with pointers, and manipulated with available
   * functions.
   *
   * Remember to include always <code>mps/mps.h</code> and not the other headers in the
   * directory. This is not supported and not guaranteed to work.
   *
   * For example you can instantiate a pointer to a new <code>mps_context</code> with
   * the function <code>mps_context_new()</code>, add the input to it, solve the polynoial
   * with <code>mps_mpsolve()</code> and then free its resources with <code>mps_context_free()</code>.
   *
   * The data structures that will be mostly used are:
   * -# <code>mps_context</code>: This structure holds the state of the computation and must
   * be instanciated for every polynomial solving operation. After allocating a pointer to
   * an <code>mps_context</code> you should, generally,:
   *   -# Set the input data and output requirements (i.e. the input polynomial, the desired
   *   output digits and goal, ...)
   *   -# Ask libmps to solve the polynomial by calling <code>mps_mpsolve()</code>
   *   -# Retrieve the computed data with the appropriate accessors functions
   * <code>mps_context_get_roots_d()</code> or <code>mps_context_get_roots_m()</code>.
   * All the functions usable on an <code>mps_context</code> pointer are available in
   * status.h.
   *
   * -# <code>mps_monomial_poly</code>: A polynomial given by its coefficients. Can be allocated
   * with <code>mps_monomial_poly_new()</code> and manipulated with the functions in 
   * monomial-poly.h. Once it is the desired polynomial to solve you can call
   * <code>mps_context_set_input_poly()</code> to set it as the active polynomial to solve.
   *
   * -# <code>mps_secular_equation</code>: The same as the monomial poly, but for secular equations.
   * See secular-equation.h for some functions to allocate, free and manipulate them.
   *
   * @subsection async Calling MPSolve asynchronously
   * 
   * Desktop application using MPSolve to solve polynomials may want to call a 
   * non-blocking version of mps_mpsolve(). This is provided inside the package
   * as mps_mpsolve_async(). 
   *
   * This routine will return a <code>mps_handle</code> pointer that can be used to wait
   * for the result by calling mps_mpsolve_wait() on it.
   */

  /*
   *    ====== ROUTINES EXPOSED TO THE INTERFACE ======
   */

  /* functions in mps_defaults.c */
  void mps_set_default_values (mps_context * s);

  /* Functions in mps_main.c */
  void mps_mpsolve (mps_context * s);
  void mps_standard_mpsolve (mps_context * s);

  /* functions in mps_interface.c */
  void * mps_malloc (size_t size);
  void * mps_realloc (void * pointer, size_t size);

  void mps_mpsolve_async (mps_context * s, mps_callback callback, void * user_data);

  /* Macros to init pointer and/or vectors in a convenient way */
#define mps_new(type) ((type *) mps_malloc (sizeof (type)))
#define mps_newv(type, n) ((type *) mps_malloc (sizeof (type) * (n)))

#ifdef  __cplusplus
}
#endif

#ifdef __UNDEF_CPLUSPLUS
}
#endif

#endif                          /* MPS_INTERFACE_H */
