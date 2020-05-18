/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/**
 * @file
 *
 * This file is the header for the libmps library. Including
 * this file is needed to access all the MPSolve routines by
 * MPSolve internals.
 *
 * @brief Header file for libmps
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef MPS_CORE_H_
#define MPS_CORE_H_

#ifdef __cplusplus
#define __MPS_NOT_DEFINE_BOOL
#endif

#ifdef __MPS_MATLAB_MODE
#define __MPS_NOT_DEFINE_BOOL
#endif

#ifdef HAVE_HIDDEN_VISIBILITY_ATTRIBUTE
  #ifdef MPS_PUBLISH_PRIVATE_METHODS
    #define MPS_PRIVATE
  #else
    #define MPS_PRIVATE __attribute__((visibility ("hidden")))
  #endif
#else
  #define MPS_PRIVATE
#endif

#ifdef __cplusplus
  #define MPS_BEGIN_DECLS extern "C" {
  #define MPS_END_DECLS }
#else
  #define MPS_BEGIN_DECLS
  #define MPS_END_DECLS
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* This header should be included first since it contains all the forward
 * declaration that must be available in the others. */
#include <mps/types.h>

/* Include more types that are needed in the declarations of the generic
 * functions, such as DPE and MP. */
#include <mps/mt.h>
#include <mps/gmptools.h>
#include <mps/mpcf.h>
#include <mps/link.h>
#include <mps/polynomial.h>

/* Public types, mostly custom polynomial types such as Chebyshev, Monomial and
 * Secular equations. */
#include <mps/matrix.h>
#include <mps/chebyshev.h>
#include <mps/monomial-matrix-poly.h>
#include <mps/monomial-poly.h>
#include <mps/secular-equation.h>
#include <mps/nroots-polynomial.h>
#include <mps/regeneration-driver.h>

/* Public interface functions for MPSolve */
#include <mps/approximation.h>
#include <mps/context.h>
#include <mps/debug.h>
#include <mps/interface.h>
#include <mps/parser.h>

/* Private inclusions. Please note that these header files may not be distributed with
 * MPSolve, so it's safe to use them only for internal functions. */
#ifdef _MPS_PRIVATE

#ifndef getline
MPS_BEGIN_DECLS
ssize_t getline (char **lineptr, size_t *n, FILE *stream);
MPS_END_DECLS
#endif
#include <mps/private/system/abstract-input-stream.h>
#include <mps/private/system/file-input-stream.h>
#include <mps/private/system/memory-file-stream.h>
#include <mps/private/aberth.h>
#include <mps/private/algorithms.h>
#include <mps/private/cluster.h>
#include <mps/private/convex.h>
#include <mps/private/data.h>
#include <mps/private/hessenberg-determinant.h>
#include <mps/private/horner.h>
#include <mps/private/jacobi-aberth.h>
#include <mps/private/improve.h>
#include <mps/private/input-buffer.h>
#include <mps/private/input-output.h>
#include <mps/private/list.h>
#include <mps/private/mandelbrot-user.h>
#include <mps/private/newton.h>
#include <mps/private/options.h>
#include <mps/private/radii.h>
#include <mps/private/secular-evaluation.h>
#include <mps/private/solve.h>
#include <mps/private/sort.h>
#include <mps/private/starting.h>
#include <mps/private/starting-configuration.h>
#include <mps/private/threading.h>
#include <mps/private/tools.h>
#include <mps/private/touch.h>
#include <mps/private/utils.h>
#include <mps/private/formal/formal-monomial.h>
#include <mps/private/formal/formal-polynomial.h>
#include <mps/private/secular-regeneration.h>
#endif

#endif                          /* ndef MPSCORE_H */
