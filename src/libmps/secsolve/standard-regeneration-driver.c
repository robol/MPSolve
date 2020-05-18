/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>

MPS_PRIVATE mps_boolean
mps_standard_regeneration_driver_update_fsecular_equation (mps_context * s, 
							   mps_polynomial * p, 
							   mps_approximation ** approximations, 
							   mps_secular_equation * sec)
{
  mps_boolean successful_regeneration = true;
  cplx_t * old_a = NULL, * old_b = NULL;
  cdpe_t * old_db = NULL;
  mpcf_t * old_mb = NULL;
  int i;

  s->mpwp = MPS_SECULAR_EQUIVALENT_FP_PRECISION;

  /* Allocate old_a and old_b */
  old_mb = mpcf_valloc (s->n);
  for (i = 0; i < s->n; i++)
    mpcf_init2 (old_mb[i], approximations[i]->wp);  

  old_a = cplx_valloc (s->n);
  old_b = cplx_valloc (s->n);
  old_db = cdpe_valloc (s->n);

  /* Copy the old coefficients, and set the new
   * b_i with the current approximationss approximations. */
  for (i = 0; i < s->n; i++)
    {
      cplx_set (old_a[i], sec->afpc[i]);
      cplx_set (old_b[i], sec->bfpc[i]);
      cdpe_set_x (old_db[i], old_b[i]);
      mpcf_set_cplx (old_mb[i], old_b[i]);
      mpcf_set_prec (sec->bmpc[i], s->mpwp);
      mpcf_set (sec->bmpc[i], approximations[i]->mvalue);
    }

  /* Regeneration */
  if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp (s, old_db, old_mb)))
    {
      for (i = 0; i < s->n; i++)
	{
	  cplx_set (sec->afpc[i], old_a[i]);
	  cplx_set (sec->bfpc[i], old_b[i]);
	}
    }
  else
    {
      mps_secular_ga_update_coefficients (s);
      for (i = 0; i < s->n; i++)
	{
	  /* We may risk that NaN or inf have been introduced because of huge
	   * coefficients computed, so let's check it and in the case of failure
	   * switch to DPE. */
	  if (cplx_check_fpe (sec->afpc[i]) || cplx_check_fpe (sec->bfpc[i]) ||
	      (cplx_mod (sec->afpc[i]) > 1.0e300) ||
	      (cplx_mod (sec->bfpc[i]) > 1.0e300))
	    {
	      successful_regeneration = false;
	      if (s->debug_level & MPS_DEBUG_REGENERATION)
		MPS_DEBUG (s, "Found floating point exception in regenerated coefficients, reusing old ones.");
	      
	      for (i = 0; i < s->n; i++)
		{
		  cplx_set (sec->afpc[i], old_a[i]);
		  cplx_set (sec->bfpc[i], old_b[i]);
		}
	      break;
	    }
	  
	  if (s->debug_level & MPS_DEBUG_REGENERATION)
	    {
	      MPS_DEBUG_CPLX (s, sec->afpc[i], "sec->afpc[%d]", i);
	      MPS_DEBUG_CPLX (s, sec->bfpc[i], "sec->bfpc[%d]", i);
	    }
	  
	  mpcf_set_cplx (approximations[i]->mvalue, approximations[i]->fvalue);
	}
    }
  
  cplx_vfree (old_a);
  cplx_vfree (old_b);
  cdpe_vfree (old_db);  

  mpcf_vclear (old_mb, s->n);
  mpcf_vfree (old_mb);  

  return successful_regeneration;
}

MPS_PRIVATE mps_boolean
mps_standard_regeneration_driver_update_dsecular_equation (mps_context * s, 
							   mps_polynomial * p, 
							   mps_approximation ** approximations, 
							   mps_secular_equation * sec)
{
  mps_boolean successful_regeneration = true;
  int i;
  cdpe_t * old_da = NULL;
  cdpe_t * old_db = NULL;
  mpcf_t * old_mb = NULL;

  old_mb = mpcf_valloc (s->n);
  for (i = 0; i < s->n; i++)
    mpcf_init2 (old_mb[i], approximations[i]->wp);  

  s->mpwp = MPS_SECULAR_EQUIVALENT_FP_PRECISION;

  /* Allocate old_a and old_b */
  old_da = cdpe_valloc (s->n);
  old_db = cdpe_valloc (s->n);

  /* Copy the old coefficients, and set the new
   * b_i with the current approximationss approximations. */
  for (i = 0; i < s->n; i++)
    {
      cdpe_set (old_da[i], sec->adpc[i]);
      cdpe_set (old_db[i], sec->bdpc[i]);
      mpcf_get_cdpe (sec->bdpc[i], approximations[i]->mvalue);
      mpcf_set_cdpe (old_mb[i], old_db[i]);
      mpcf_set_prec (sec->bmpc[i], s->mpwp);
      mpcf_set (sec->bmpc[i], approximations[i]->mvalue);
    }

  /* Regeneration */
  if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp (s, old_db, old_mb)))
    {
      MPS_DEBUG (s, "Regeneration failed");
      for (i = 0; i < s->n; i++)
	{
	  cdpe_set (sec->adpc[i], old_da[i]);
	  cdpe_set (sec->bdpc[i], old_db[i]);
	  mpcf_set_cdpe (old_mb[i], old_db[i]);
	  mpcf_set_cdpe (sec->ampc[i], old_da[i]);
	  mpcf_set_cdpe (sec->bmpc[i], old_db[i]);
	}

      mps_secular_ga_update_coefficients (s);
      successful_regeneration = false;
      goto cleanup;
    }
  else
    {
      mps_secular_ga_update_coefficients (s);
      for (i = 0; i < s->n; i++)
	mpcf_set_cdpe (approximations[i]->mvalue, approximations[i]->dvalue);
      mps_secular_set_radii (s);
    }

  if (s->debug_level & MPS_DEBUG_REGENERATION)
    {
      for (i = 0; i < s->n; i++)
	{
	  MPS_DEBUG_CDPE (s, sec->bdpc[i], "sec->bdpc[%d]", i);
	  MPS_DEBUG_CDPE (s, sec->adpc[i], "sec->adpc[%d]", i);
	}
    }

  /* Free data */
 cleanup:

  cdpe_vfree (old_da);
  cdpe_vfree (old_db);

  mpcf_vclear (old_mb, MPS_POLYNOMIAL (sec)->degree);
  mpcf_vfree  (old_mb);
      
  return successful_regeneration;
}

MPS_PRIVATE mps_boolean
mps_standard_regeneration_driver_update_msecular_equation (mps_context * s, 
						       mps_polynomial * p,
						       mps_approximation ** approximations, 
						       mps_secular_equation * sec)
{
  mps_boolean successful_regeneration = true;
  int i;
  mpcf_t* old_ma = NULL, *old_mb = NULL;
  cdpe_t * old_db = NULL;

  /* Allocate old_a and old_b */
  old_ma = mpcf_valloc (s->n);
  old_mb = mpcf_valloc (s->n);
  old_db = cdpe_valloc (s->n);

  mpcf_vinit2 (old_ma, s->n, s->mpwp);
  mpcf_vinit2 (old_mb, s->n, s->mpwp);

  /* Copy the old coefficients, and set the new
   * b_i with the current roots approximations. */
  for (i = 0; i < s->n; i++)
    {
      mpcf_set (old_ma[i], sec->ampc[i]);
      mpcf_set (old_mb[i], sec->bmpc[i]);
      mpcf_set_prec (sec->bmpc[i], mpcf_get_prec (s->root[i]->mvalue));
      mpcf_set (sec->bmpc[i], s->root[i]->mvalue);
      mpcf_get_cdpe (old_db[i], old_mb[i]);
    }

  /* Regeneration */
  if ((successful_regeneration = mps_secular_ga_regenerate_coefficients_mp (s, old_db, old_mb)))
    {
      mps_secular_ga_update_coefficients (s);
      /* Finally set radius according to new computed a_i coefficients,
       * if they are convenient   */
      mps_secular_set_radii (s);
    }
  else
    MPS_DEBUG (s, "Regeneration failed");

  if (s->debug_level & MPS_DEBUG_REGENERATION)
    {
      MPS_DEBUG (s, "Dumping regenerated coefficients");
      for (i = 0; i < s->n; i++)
	{
	  MPS_DEBUG_MPC (s, s->mpwp, sec->ampc[i], "ampc[%d]", i);
	  MPS_DEBUG_MPC (s, s->mpwp, sec->bmpc[i], "bmpc[%d]", i);
	}
    }

  mpcf_vclear (old_ma, s->n);
  mpcf_vfree (old_ma);
  rdpe_vfree (old_db);

  return successful_regeneration;
}

static mps_regeneration_driver _mps_standard_regeneration_driver_instance = { 
  mps_standard_regeneration_driver_update_fsecular_equation, 
  mps_standard_regeneration_driver_update_dsecular_equation,
  mps_standard_regeneration_driver_update_msecular_equation,
  NULL
};

mps_regeneration_driver *
mps_regeneration_driver_new_standard (mps_context * ctx)
{
  return &_mps_standard_regeneration_driver_instance;
}

void
mps_regeneration_driver_free (mps_context * ctx, mps_regeneration_driver * rd)
{
  if (rd->free)
    rd->free(ctx, rd);

  free (rd);
}
