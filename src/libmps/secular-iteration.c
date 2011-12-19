#include <mps/mps.h>
#include <math.h>
#include <string.h>

/* Routine called to perform the floating point newton iterations. */
void *
__mps_secular_ga_fiterate_worker (void* data_ptr)
{
  mps_thread_worker_data *data = (mps_thread_worker_data *) data_ptr;
  mps_status *s = data->s;
  /* mps_secular_equation *sec = s->secular_equation; */
  int i;
  cplx_t corr, abcorr, froot;
  double modcorr;
  mps_thread_job job;

  mps_secular_iteration_data it_data;

  while ((*data->nzeros < s->n))
    {
      job = mps_thread_job_queue_next (s, data->queue);
      i = job.i;

      if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
          return NULL;
        }

      if (s->again[i])
	{
          /* Lock this roots to make sure that we are the only one working on it */
          pthread_mutex_lock (&data->roots_mutex[i]);
	  cplx_set (froot, s->froot[i]);

	  (*data->it)++;

	  it_data.k = i;
	  mps_secular_fnewton (s, froot, &s->frad[i], corr,
			       &s->again[i], &it_data, false);

	  /* Apply Aberth correction */
	  mps_faberth (s, i, abcorr);
	  cplx_mul_eq (abcorr, corr);
	  cplx_sub (abcorr, cplx_one, abcorr);
	  cplx_div (abcorr, corr, abcorr);
	  cplx_sub_eq (froot, abcorr);

	  /* Correct the radius */
	  modcorr = cplx_mod (abcorr);
	  s->frad[i] += modcorr;

	  if (modcorr < cplx_mod (froot) * 8.0 * DBL_EPSILON)
	    {
	      if (s->debug_level & MPS_DEBUG_PACKETS)
		{
		  MPS_DEBUG (s, "Setting again to false on root %d because aberth correction is less than precision", i);
		}
	      s->again[i] = false;
	    }

	  if (!s->again[i])
	    {
	      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
		MPS_DEBUG (s, "Root %d again was set to false on iteration %d", i, *data->it);
	      (*data->nzeros)++;
	    }

	  cplx_set (s->froot[i], froot);
	  pthread_mutex_unlock (&data->roots_mutex[i]);
	}
    }

  return NULL;
}

/**
 * @brief Routine that performs a block of iteration
 * in floating point on the secular equation.
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_fiterate (mps_status * s, int maxit, mps_boolean just_regenerated)
{
  int computed_roots = 0;
  int i;
  int nit = 0;
  int it_threshold;
  double * fradii = double_valloc (s->n);

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  mps_secular_equation *sec = s->secular_equation;

  mps_thread_worker_data *data;
  pthread_mutex_t *aberth_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);
  pthread_mutex_t *roots_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);

  for (i = 0; i < s->n; i++)
    {
      pthread_mutex_init (roots_mutex + i, NULL);
      pthread_mutex_init (aberth_mutex + i, NULL);
    }

  data = mps_newv (mps_thread_worker_data, s->n_threads);

  MPS_DEBUG_THIS_CALL;

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

  mps_thread_job_queue *queue = mps_thread_job_queue_new (s);
  queue->max_iter = it_threshold;

  for (i = 0; i < s->n_threads; i++)
    {
      data[i].it = &nit;
      data[i].nzeros = &computed_roots;
      data[i].s = s;
      data[i].thread = i;
      data[i].n_threads = s->n_threads;
      data[i].aberth_mutex = aberth_mutex;
      data[i].roots_mutex = roots_mutex;
      data[i].queue = queue;

      mps_thread_pool_assign (s, s->pool, __mps_secular_ga_fiterate_worker, data + i);
    }

  mps_thread_pool_wait (s, s->pool);

  /* Check if the roots are improvable in floating point */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations with floating point arithmetic",
                       nit);

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
      mps_dump (s);

  if (nit <= it_threshold && just_regenerated)
    {
      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	{
	  MPS_DEBUG (s, "Setting approximation as best_approx");
	}
      s->secular_equation->best_approx = true;
    }

  /* Compute the inclusion radii with Gerschgorin so we can compute
   * clusterizations for the roots. */
  mps_fradii (s, fradii);
  mps_fcluster (s, fradii, 2.0 * s->n); 
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

  free (fradii);

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
mps_secular_ga_diterate (mps_status * s, int maxit, mps_boolean just_regenerated)
{
  int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;
  int it_threshold;
  mps_secular_equation *sec = s->secular_equation;
  mps_secular_iteration_data data;
  rdpe_t * drad = rdpe_valloc (s->n);
  rdpe_t * old_radii = rdpe_valloc (s->n);

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  MPS_DEBUG_THIS_CALL;

  sec->best_approx = false;

  memcpy (old_radii, s->drad, s->n * sizeof (__rdpe_struct));

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
  /* mps_dump (s); */
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
	      rdpe_t rtmp;

	      if (cdpe_eq (s->droot[i], sec->bdpc[i]))
		continue;
              nit++;

	      /* Prepare data for the dnewton routine */
	      data.k = i;

              mps_secular_dnewton (s, s->droot[i], s->drad[i], corr,
                                   &s->again[i], &data, false);
	      
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

	      cdpe_mod (rtmp, s->droot[i]);
	      rdpe_mul_eq_d (rtmp, 8.0 * DBL_EPSILON);
	      if (rdpe_lt (modcorr, rtmp))
		{
		  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
		    MPS_DEBUG (s, "Setting again to false on root %d because Aberth correction is less than machine epsilon", i);
		  s->again[i] = false;
		  
		}

              if (!s->again[i])
		{
		  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
		    MPS_DEBUG (s, "Root %d again was set to false on iteration %d", i, iterations);
		  computed_roots++;
		}
            }
        }
    }

  /* Check if no more than 2 iterations per root
   * were computed, and in that case state that
   * a coefficient regeneration won't be of much help */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations", nit);
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      mps_dump (s);
    }

  if (nit <= it_threshold && just_regenerated)
    {
      if (s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting best_approx to true");
      s->secular_equation->best_approx = true;
    }

  mps_dradii (s, drad);
  mps_dcluster (s, drad, 2.0 * s->n);
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

  rdpe_vfree (old_radii);
  rdpe_vfree (drad);

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
mps_secular_ga_miterate (mps_status * s, int maxit, mps_boolean just_regenerated)
{
  int computed_roots = 0;
  int iterations = 0;
  int i, k;
  int nit = 0;
  int it_threshold;

  mpc_t corr, abcorr;
  cdpe_t ctmp;
  rdpe_t modcorr, rtmp;
  rdpe_t * drad = rdpe_valloc (s->n);

  mps_secular_equation * sec = s->secular_equation;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  /* The data used to determined if the radius has been
   * set and to intercommunicate with the iterator */
  mps_secular_iteration_data user_data;

  /* Cluster data used in iterations */
  mps_cluster_item * c_item;
  mps_cluster * cluster;
  mps_root * root;

  MPS_DEBUG_WITH_INFO (s, "Precision is at %ld bits", s->mpwp);
  MPS_DEBUG_RDPE (s, s->mp_epsilon, "Machine epsilon is s->mp_epsilon");

  MPS_DEBUG_THIS_CALL;

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
      MPS_DEBUG_WITH_INFO (s, "%d roots are already approximated on the start of miterate", computed_roots);
    }

  /* Set the iteration threshold to two times the remaining roots
   * to compute. */
  it_threshold = 2 * (s->n - computed_roots);

  while (computed_roots < s->n && iterations < maxit)
    {
      /* Increase iterations counter */
      iterations++;

      for (c_item = s->clusterization->first; c_item != NULL; c_item = c_item->next)
        {
	  cluster = c_item->cluster;
          for (root = cluster->first; root != NULL; root = root->next)
            {
              k = root->k;
              if (s->again[k])
                {
                  nit++;

		  /* Set the correct index */
		  user_data.k = k;
                  mps_secular_mnewton (s, s->mroot[k], s->drad[k], corr,
                                       &s->again[k], &user_data, false);

                  /* Apply Aberth correction */
                  mps_maberth_s (s, k, cluster, abcorr);
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

  if (nit <= it_threshold && just_regenerated)
    s->secular_equation->best_approx = true;

  /* Perform cluster analysis */
  mps_mradii (s, drad);
  mps_mcluster (s, drad, 2.0 * s->n);
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
    }
  
  /* Clock the routine */
#ifndef DISABLE_DEBUG
  s->mp_iteration_time += mps_stop_timer (my_clock);
#endif

  rdpe_vfree (drad);

  /* Return the number of approximated roots */
  return computed_roots;
}
