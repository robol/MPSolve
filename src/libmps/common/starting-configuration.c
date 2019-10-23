/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#include <mps/mps.h>

void 
mps_starting_configuration_clear (mps_context * ctx, mps_starting_configuration * c)
{
  if (c->fradii)
    free (c->fradii);

  if (c->dradii)
    free (c->dradii);

  if (c->partitioning)
    free (c->partitioning); 

  c->fradii = NULL; 
  c->dradii = NULL; 
  c->partitioning = NULL;

  c->n_radii = 0;
}
