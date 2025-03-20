/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#include <mps/mps.h>
#include <string.h>
#include <ctype.h>

/**
 * @brief Return a string representation of MPSolve's version.
 *
 * @return A pointer to a const string containing the version number in MPSolve.
 */
const char *
mps_get_version ()
{
  return PACKAGE_VERSION;
}

/**
 * @brief Return a string representation of MPSolve's major version.
 *
 * @return An unsigned integer with the major version.
 */
unsigned int
mps_get_major_version()
{
  return MPS_MAJOR_VERSION;
}

/**
 * @brief Return a string representation of MPSolve's minor version.
 *
 * @return An unsigned integer with the minor version.
 */
unsigned int
mps_get_minor_version()
{
  return MPS_MINOR_VERSION;
}

/**
 * @brief Return a string representation of MPSolve's patch version.
 *
 * @return An unsigned integer with the patch version.
 */
unsigned int
mps_get_patch_version()
{
  return MPS_PATCH_VERSION;
}
