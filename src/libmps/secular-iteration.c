#include <mps/debug.h>
#include <mps/core.h>
#include <mps/link.h>
#include <mps/secular.h>
#include <mps/debug.h>
#include <math.h>
#include <string.h>


/**
 * @brief Routine that performs a block of iteration
 * in floating point on the secular equation.
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_fiterate (mps_status * s, int maxit)
{
  MPS_DEBUG_THIS_CALL;

  int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;
  int it_threshold;
  mps_secular_iteration_data data;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  mps_secular_equation *sec = s->secular_equation;

  double old_rad;
  cplx_t old_root;

  sec->best_approx = false;

  /* Mark the approximated roots as ready for output */
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated */
      if (s->status[i][0] == 'a' || s->status[i][0] == 'i'
	  || s->status[i][0] == 'o')
	{
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    {
	      MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
	    }
	  s->again[i] = false;
	}

      if (!s->again[i])
        computed_roots++;
    }

  /* Set the iterations threshold to 2 iterations
   * for every non approximated root. */
  it_threshold = 2 * (s->n - computed_roots);

  if (s->debug_level & MPS_DEBUG_PACKETS)
    {
      MPS_DEBUG (s, "There are %d roots with again set to false", computed_roots);
    }

  while (computed_roots < s->n && iterations < maxit)
    {
      cplx_t corr, abcorr;
      double modcorr;

      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->n; i++)
        {
          if (s->again[i])
            {
	      /* MPS_DEBUG (s, "Iterating on root %d", i); */
	      /* if (cplx_eq (s->froot[i], sec->bfpc[i])) */
	      /* 	continue; */

              nit++;

	      /* Save the old root in case that we occur in floating point
	       * exceptions */
              cplx_set (old_root, s->froot[i]);
              old_rad = s->frad[i];
	      
	      /* Prepare the data to be passed for secular-newton */
	      data.k = i;

              mps_secular_fnewton (s, s->froot[i], &s->frad[i], corr,
                                   &s->again[i], &data);

              /* Apply Aberth correction */
              mps_faberth (s, i, abcorr);
              cplx_mul_eq (abcorr, corr);
              cplx_sub (abcorr, cplx_one, abcorr);
              cplx_div (abcorr, corr, abcorr);
              cplx_sub_eq (s->froot[i], abcorr);

              /* Check if we need to switch to DPE */
              if (isnan (cplx_Re (s->froot[i]))
                  || isinf (cplx_Re (s->froot[i]))
                  || isnan (cplx_Im (s->froot[i]))
                  || isinf (cplx_Im (s->froot[i])) || isnan (s->frad[i])
                  || isinf (s->frad[i])
		  || s->status[i][0] == 'x')
                {
		  if (s->status[i][0] != 'x')
		    {
		      MPS_DEBUG_WITH_INFO (s,
					   "Switching to DPE phase because NAN or INF was introduced in computation");
		    }
		  else
		    {
		      s->status[i][0] = 'c';
		      MPS_DEBUG_WITH_INFO (s, "Switching to DPE phase because there is an approximation not representable in double");
		    }
                  cplx_set (s->froot[i], old_root);
                  s->frad[i] = old_rad;
                  s->lastphase = dpe_phase;

                  /* Copy roots, radius and coefficients */
                  for (i = 0; i < s->n; i++)
                    {
                      cdpe_set_x (s->droot[i], s->froot[i]);
                      rdpe_set_d (s->drad[i], s->frad[i]);

		      cdpe_set_x (sec->adpc[i], sec->afpc[i]);
		      cdpe_set_x (sec->bdpc[i], sec->bfpc[i]);

		      MPS_DEBUG_CDPE (s, sec->adpc[i], "sec->adpc[%d]", i);
		      MPS_DEBUG_CDPE (s, sec->bdpc[i], "sec->bdpc[%d]", i);
                    }

#ifndef DISABLE_DEBUG
                  s->fp_iteration_time += mps_stop_timer (my_clock);
#endif
                  return -1;
                }

              /* Correct the radius */
              modcorr = cplx_mod (abcorr);
              s->frad[i] += modcorr;

	      /* MPS_DEBUG (s, "modcorr = %e", modcorr); */
	      /* MPS_DEBUG_CPLX (s, s->froot[i], "s->froot[%d]", i); */

	      if (modcorr < cplx_mod (s->froot[i]) * 4 * DBL_EPSILON)
		{
		  if (s->debug_level & MPS_DEBUG_PACKETS)
		    {
		      MPS_DEBUG (s, "Setting again to false on root %d because aberth correction is less than precision", i);
		    }
		  s->again[i] = false;
		}

              if (!s->again[i])
                computed_roots++;
            }
        }
    }

  /* Check if the roots are improvable in floating point */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations with floating point arithmetic",
                       nit);

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
      mps_dump (s, s->logstr);

  if (nit <= it_threshold)
    {
      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	{
	  MPS_DEBUG (s, "Setting approximation as best_approx");
	}
      s->secular_equation->best_approx = true;
    }

   mps_fcluster (s, 2.0 * s->n); 
   mps_fmodify (s, false); 

  /* These lines are used to debug the again vector, but are not useful
   * at the moment being */
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      __MPS_DEBUG (s, "Again vector = ");
      for(i = 0; i < s->n; i++)
	{
	  fprintf (s->logstr, "%d ", s->again[i]);
	}
      fprintf (s->logstr, "\n");
    }

  /* Count time taken  */
#ifndef DISABLE_DEBUG
  s->fp_iteration_time += mps_stop_timer (my_clock);
#endif

  /* Return the number of approximated roots */
  return computed_roots;
}

/**
 * @brief Routine that performs a block of iteration
 * in DPE on the secular equation.
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_diterate (mps_status * s, int maxit)
{
  MPS_DEBUG_THIS_CALL;

  int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;
  int it_threshold;
  mps_secular_equation *sec = s->secular_equation;
  mps_secular_iteration_data data;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  sec->best_approx = false;

  /* Iterate with newton until we have good approximations
   * of the roots */
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated */
      if (s->status[i][0] == 'a' || s->status[i][0] == 'i'
	  || s->status[i][0] == 'o')
	{
	  MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
	  s->again[i] = false;
	}
    }

  for (i = 0; i < s->n; i++)
    {
      if (!s->again[i])
        computed_roots++;
    }

  /* Set the iterations threshold to 2 iterations
   * for every non approximated root. */
  it_threshold = 2 * (s->n - computed_roots);

  if (s->debug_level & MPS_DEBUG_PACKETS)
    MPS_DEBUG (s, "There are %d roots with again set to false", computed_roots);

  /* Use this dump only for debugging purpose */
  /* mps_dump (s, s->logstr); */
  while (computed_roots < s->n && iterations < maxit)
    {
      cdpe_t corr, abcorr;
      rdpe_t modcorr;

      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->n; i++)
        {
          if (s->again[i])
            {
	      if (cdpe_eq (s->droot[i], sec->bdpc[i]))
		continue;
              nit++;

	      /* Prepare data for the dnewton routine */
	      data.k = i;

              mps_secular_dnewton (s, s->droot[i], s->drad[i], corr,
                                   &s->again[i], &data);
	      
              /* Apply Aberth correction */
              mps_daberth (s, i, abcorr);
              cdpe_mul_eq (abcorr, corr);
              cdpe_sub (abcorr, cdpe_one, abcorr);

              if (cdpe_ne (abcorr, cdpe_zero))
                {
                  cdpe_div (abcorr, corr, abcorr);
                  cdpe_sub_eq (s->droot[i], abcorr);

                  /* Correct the radius */
                  cdpe_mod (modcorr, abcorr);
                  rdpe_add_eq (s->drad[i], modcorr);
                }
              else
                  s->again[i] = true;

              if (!s->again[i])
                computed_roots++;
            }
        }
    }

  /* Check if no more than 2 iterations per root
   * were computed, and in that case state that
   * a coefficient regeneration won't be of much help */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations", nit);
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      mps_dump (s, s->logstr);
    }

  if (nit <= it_threshold)
    {
      if (s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting best_approx to true");
      s->secular_equation->best_approx = true;
    }

  mps_dcluster (s, 2.0 * s->n);
  mps_dmodify (s, false);

  /* These lines are used to debug the again vector, but are not useful
   * at the moment being */
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      __MPS_DEBUG (s, "Again vector = ");
      for(i = 0; i < s->n; i++)
	{
	  fprintf (s->logstr, "%d ", s->again[i]);
	}
      fprintf (s->logstr, "\n");
    }

  /* Clock the routine */
#ifndef DISABLE_DEBUG
  s->dpe_iteration_time += mps_stop_timer (my_clock);
#endif

  /* Return the number of approximated roots */
  return computed_roots;
}

/**
 * @brief Routine that performs a block of iteration
 * in Multiprecision on the secular equation.
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_miterate (mps_status * s, int maxit)
{
  MPS_DEBUG_THIS_CALL;

  MPS_DEBUG_WITH_INFO (s, "Precision is at %ld bits", s->mpwp);
  MPS_DEBUG_RDPE (s, s->mp_epsilon, "Machine epsilon is s->mp_epsilon");

  int computed_roots = 0;
  int iterations = 0;
  int i, j, k;
  int nit = 0;
  int it_threshold;
  int old_cr;

  mpc_t corr, abcorr;
  cdpe_t ctmp;
  rdpe_t modcorr, rtmp;

  mps_secular_equation * sec = s->secular_equation;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  /* The data used to determined if the radius has been
   * set and to intercommunicate with the iterator */
  mps_secular_iteration_data user_data;

  /* Init data with the right precision */
  mpc_init2 (corr, s->mpwp);
  mpc_init2 (abcorr, s->mpwp);

  sec->best_approx = false;

  /* Iterate with newton until we have good approximations
   * of the roots */
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated */
      if (s->status[i][0] == 'a' || s->status[i][0] == 'i'
	  || s->status[i][0] == 'o')
	{
	  MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
	  s->again[i] = false;
	}

      if (s->again[i])
        s->rootwp[i] = s->mpwp;
      else
        computed_roots++;
    }

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS) 
    {
      MPS_DEBUG_WITH_INFO (s, "%d roots are already approximated on the start of miterate", computed_roots)
    }

  /* Set the iteration threshold to two times the remaining roots
   * to compute. */
  old_cr = computed_roots;
  it_threshold = 2 * (s->n - computed_roots);

  while (computed_roots < s->n && iterations < maxit)
    {
      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->nclust; i++)
        {
          for (j = s->punt[i]; j < s->punt[i + 1]; j++)
            {
              k = s->clust[j];
              if (s->again[k])
                {
                  nit++;

		  /* Set the correct index */
		  user_data.k = k;
                  mps_secular_mnewton (s, s->mroot[k], s->drad[k], corr,
                                       &s->again[k], &user_data);

                  /* Apply Aberth correction */
                  mps_maberth_s (s, k, i, abcorr);
                  mpc_mul_eq (abcorr, corr);
                  mpc_ui_sub (abcorr, 1, 0, abcorr);
                  mpc_div (abcorr, corr, abcorr);
                  mpc_sub_eq (s->mroot[k], abcorr);

                  /* Correct the radius */
                  mpc_get_cdpe (ctmp, abcorr);
                  cdpe_mod (modcorr, ctmp);
                  rdpe_add_eq (s->drad[k], modcorr);

                  mpc_get_cdpe (ctmp, s->mroot[k]);
                  cdpe_mod (rtmp, ctmp);
                  rdpe_div_eq (modcorr, rtmp);

                  if (!s->again[k])
                    computed_roots++;
                }
            }
        }
    }

  /* Deallocate multiprecision local variables */
  mpc_clear (abcorr);
  mpc_clear (corr);

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      MPS_DEBUG (s, "Performed %d iterations", nit);
    }

  if (nit <= it_threshold)
    s->secular_equation->best_approx = true;

  /* Perform cluster analysis */
  mps_mcluster (s, 2.0 * s->n);
  mps_mmodify (s, false);

  /* These lines are used to debug the again vector, but are not useful
   * at the moment being */
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      __MPS_DEBUG (s, "Again vector = ");
      for (i = 0; i < s->n; i++)
	{
	  fprintf (s->logstr, "%d ", s->again[i]);
	}
      fprintf (s->logstr, "\n");
      __MPS_DEBUG (s, "Status = ");
      for (i = 0; i < s->n; i++)
	{
	  fprintf (s->logstr, "%c ", s->status[i][0]);
	}
      fprintf (s->logstr, "\n");
      // mps_dump (s, s->logstr);
    }
  
  /* Clock the routine */
#ifndef DISABLE_DEBUG
  s->mp_iteration_time += mps_stop_timer (my_clock);
#endif

  /* Return the number of approximated roots */
  return computed_roots;
}