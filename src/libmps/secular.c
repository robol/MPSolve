/*
 * mps_secular.c
 *
 *  Created on: 10/apr/2011
 *      Author: leonardo
 */

#include <mps/core.h>
#include <mps/secular.h>
#include <mps/debug.h>
#include <mps/mt.h>
#include <float.h>
#include <mps/mpc.h>


/**
 * @brief Deflate a secular equation lowering the degree of the
 * polynomial that represent it, if that is possible.
 */
void
mps_secular_deflate(mps_secular_equation* sec)
{
  int i, j, k;
  for (i = 0; i < sec->n; i++)
    {
      for (j = i + 1; j < sec->n; j++)
        {
          if (cplx_eq(sec->bfpc[i], sec->bfpc[j]))
            {
              cplx_add_eq(sec->afpc[i], sec->afpc[j]);

              /* Copy other coefficients back of one position */
              for (k = j; k < sec->n - 1; k++)
                {
                  cplx_set(sec->afpc[j], sec->afpc[j + 1]);
                  cplx_set(sec->bfpc[j], sec->bfpc[j + 1]);
                }

              /* Decrement number of coefficients */
              sec->n--;
            }
        }
    }
}

/**
 * @brief Create a new secular equation struct
 */
mps_secular_equation*
mps_secular_equation_new(cplx_t* afpc, cplx_t* bfpc, unsigned long int n)
{

  int i;

  /* Allocate the space for the new struct */
  mps_secular_equation* s = (mps_secular_equation*) malloc(
      sizeof(mps_secular_equation));

  /* Copy data in the struct, so the user shall not worry about the scope of
   * its input data. */
  s->afpc = cplx_valloc(n);
  s->bfpc = cplx_valloc(n);

  /* Allocate complex dpe coefficients of the secular equation */
  s->adpc = cdpe_valloc(n);
  s->bdpc = cdpe_valloc(n);

  /* Allocate multiprecision complex coefficients of the secular equation */
  s->ampc = mpc_valloc(n);
  s->bmpc = mpc_valloc(n);

  /* Copy the complex coefficients passed as argument */
  for (i = 0; i < n; i++)
    {
      /* a_i coefficients */
      cplx_set(s->afpc[i], afpc[i]);

      /* b_i coefficients */
      cplx_set(s->bfpc[i], bfpc[i]);
    }

  s->n = n;
  mps_secular_deflate(s);

  for (i = 0; i < s->n; i++)
    {
      cdpe_init(s->adpc[i]);
      cdpe_set_x(s->adpc[i], s->afpc[i]);

      mpc_init(s->ampc[i]);
      mpc_set_cplx(s->ampc[i], s->afpc[i]);

      cdpe_init(s->bdpc[i]);
      cdpe_set_x(s->bdpc[i], s->bfpc[i]);

      mpc_init(s->bmpc[i]);
      mpc_set_cplx(s->bmpc[i], s->bfpc[i]);
    }

  return s;
}

void
mps_secular_equation_free(mps_secular_equation* s)
{
  /* Free internal data */
  cplx_vfree(s->afpc);
  cplx_vfree(s->bfpc);

  /* ...and then release it */
  free(s);
}



/**
 * @brief Evaluate secular equation in the point x.
 */
void
mps_secular_evaluate(mps_status* s, cplx_t x, cplx_t sec_ev)
{
  cplx_t ctmp;
  int i;
  mps_secular_equation* sec = (mps_secular_equation*) s->user_data;
  cplx_set(sec_ev, cplx_zero);

  for (i = 0; i < s->n; i++)
    {
      /* Compute 1 / (x - b_i) */
      cplx_sub(ctmp, x, sec->bfpc[i]);
      cplx_inv_eq(ctmp);

      /* Compute a_i / (x - b_i) */
      cplx_mul_eq(ctmp, sec->afpc[i]);

      /* Sum to the secular eqation */
      cplx_add_eq(sec_ev, ctmp);
    }

  cplx_sub_eq(sec_ev, cplx_one);
}


void
mps_secular_check_data(mps_status* s, char* which_case)
{
  /* While we can't found a good criterion to check
   * the possibility to start in pure floating point we
   * use the DPE version. */
  mps_secular_equation* sec = (mps_secular_equation*) s->user_data;
  *which_case = sec->starting_case;
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
mps_secular_ga_fiterate(mps_status* s, int maxit)
{
  MPS_DEBUG_THIS_CALL

  int computed_roots = 0;
  int iterations = 0;
  int i;

  /* Iterate with newton until we have good approximations
   * of the roots */
  /* Set again to true */
  for (i = 0; i < s->n; i++)
    s->again[i] = true;

  while (computed_roots < s->n && iterations < maxit - 1)
    {
      cplx_t corr, abcorr;
      double modcorr;

      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->n; i++)
        {
          if (s->again[i])
            {
              mps_secular_fnewton(s, s->froot[i], &s->frad[i], corr,
                  &s->again[i]);

              /* Apply Aberth correction */
              mps_faberth(s, i, abcorr);
              cplx_mul_eq(abcorr, corr);
              cplx_sub(abcorr, cplx_one, abcorr);
              cplx_div(abcorr, corr, abcorr);
              cplx_sub_eq(s->froot[i], abcorr);

              /* Correct the radius */
              modcorr = cplx_mod(abcorr);
              s->frad[i] += modcorr;

              if (!s->again[i])
                computed_roots++;
            }
        }
    }

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
mps_secular_ga_diterate(mps_status* s, int maxit)
{
  MPS_DEBUG_THIS_CALL

  int computed_roots = 0;
  int iterations = 0;
  int i;

  /* Iterate with newton until we have good approximations
   * of the roots */
  /* Set again to true */
  for (i = 0; i < s->n; i++)
    s->again[i] = true;

  while (computed_roots < s->n && iterations < maxit - 1)
    {
      cdpe_t corr, abcorr;
      rdpe_t modcorr;

      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->n; i++)
        {
          if (s->again[i])
            {
              mps_secular_dnewton(s, s->droot[i], s->drad[i], corr,
                  &s->again[i]);

              /* Apply Aberth correction */
              mps_daberth(s, i, abcorr);
              cdpe_mul_eq(abcorr, corr);
              cdpe_sub(abcorr, cdpe_one, abcorr);
              cdpe_div(abcorr, corr, abcorr);
              cdpe_sub_eq(s->droot[i], abcorr);

              /* Correct the radius */
              cdpe_mod(modcorr, abcorr);
              rdpe_add_eq(s->drad[i], modcorr);

              if (!s->again[i])
                computed_roots++;
            }
        }
    }

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
mps_secular_ga_miterate(mps_status* s, int maxit)
{
  MPS_DEBUG_THIS_CALL

  int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;

  mpc_t corr, abcorr;
  cdpe_t ctmp;
  rdpe_t modcorr, drad;

  rdpe_set_2dl(drad, 1.0, -s->prec_out);

  /* Init data with the right precision */
  mpc_init2(corr, s->mpwp);
  mpc_init2(abcorr, s->mpwp);

  /* Iterate with newton until we have good approximations
   * of the roots */
  /* Allocate again and set all to true */
  for (i = 0; i < s->n; i++)
    {
      if (rdpe_gt(s->drad[i], drad))
        s->again[i] = true;
      else
        {
          s->again[i] = false;
          computed_roots++;
        }
    }

  while (computed_roots < s->n && iterations < maxit)
    {
      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->n; i++)
        {
          if (s->again[i])
            {
              nit++;
              mps_secular_mnewton(s, s->mroot[i], s->drad[i], corr,
                  &s->again[i]);

              /* Apply Aberth correction */
              mps_maberth(s, i, abcorr);
              mpc_mul_eq(abcorr, corr);
              mpc_ui_sub(abcorr, 1, 0, abcorr);
              mpc_div(abcorr, corr, abcorr);
              mpc_sub_eq(s->mroot[i], abcorr);

              /* Correct the radius */
              mpc_get_cdpe(ctmp, abcorr);
              cdpe_mod(modcorr, ctmp);
              rdpe_add_eq(s->drad[i], modcorr);

              if (!s->again[i])
                computed_roots++;
            }
        }
    }

  /* Deallocate multiprecision local variables */
  mpc_clear(abcorr);
  mpc_clear(corr);

  MPS_DEBUG(s, "Performed %d iterations", nit)

  /* Return the number of approximated roots */
  return computed_roots;
}

/**
 * @brief Regenerate \f$a_i\f$ and \f$b_i\f$ setting
 * \f$b_i = z_i\f$, i.e. the current root approximation
 * and recomputing \f$a_i\f$ accordingly.
 */
void
mps_secular_ga_regenerate_coefficients(mps_status* s)
{
  MPS_DEBUG_THIS_CALL

  cplx_t *old_b, *old_a;
  cdpe_t *old_db, *old_da;
  mpc_t *old_ma, *old_mb;
  mps_secular_equation *sec;
  int i, j;

  /* Declaration and initialization of the multprecision
   * variables that are used only in that case */
  mpc_t prod_old_b, prod_b, sec_ev;
  mpc_t ctmp, btmp;

  mpc_init2(prod_old_b, s->mpwp);
  mpc_init2(prod_b, s->mpwp);
  mpc_init2(sec_ev, s->mpwp);
  mpc_init2(ctmp, s->mpwp);
  mpc_init2(btmp, s->mpwp);

  sec = (mps_secular_equation*) s->user_data;

  switch (s->lastphase)
    {

  /* If we are in the float phase regenerate coefficients
   * starting from floating point */
  case float_phase:

    /* Allocate old_a and old_b */
    old_a = cplx_valloc(s->n);
    old_b = cplx_valloc(s->n);

    /* Copy the old coefficients, and set the new
     * b_i with the current roots approximations. */
    for (i = 0; i < s->n; i++)
      {
        cplx_set(old_a[i], sec->afpc[i]);
        cplx_set(old_b[i], sec->bfpc[i]);
        cplx_set(sec->bfpc[i], s->froot[i]);
      }

    /* Compute the new a_i */
    for (i = 0; i < s->n; i++)
      {
        cplx_t prod_old_b, prod_b, sec_ev;
        cplx_t ctmp, btmp;
        cplx_set(prod_old_b, cplx_one);
        cplx_set(prod_b, cplx_one);
        cplx_set(sec_ev, cplx_zero);

        for (j = 0; j < sec->n; j++)
          {
            /* Compute 1 / (b - old_b) */
            cplx_sub(btmp, sec->bfpc[i], old_b[j]);

            if (cplx_eq_zero(btmp))
              {
                for (i = 0; i < sec->n; i++)
                  {
                    cplx_set(sec->afpc[i], old_a[i]);
                    cplx_set(sec->bfpc[i], old_b[i]);

                    MPS_DEBUG(s, "Cannot regenerate coefficients, reusing old ones")
                    return;
                  }
              }

            cplx_inv(ctmp, btmp);

            /* Add a_j / (b_i - old_b_j) to sec_ev */
            cplx_mul_eq(ctmp, old_a[j]);
            cplx_add_eq(sec_ev, ctmp);

            /* Multiply prod_b for
             * b_i - b_j if i \neq j and prod_old_b
             * for b_i - old_b_i.  */
            cplx_mul_eq(prod_b, btmp);
            if (i != j)
              {
                cplx_sub(ctmp, sec->bfpc[i], sec->bfpc[j]);
                cplx_mul_eq(prod_old_b, ctmp);
              }
          }

        /* Compute the new a_i as sec_ev * prod_old_b / prod_b */
        cplx_sub_eq(sec_ev, cplx_one);
        cplx_mul(sec->afpc[i], sec_ev, prod_old_b);
        cplx_div_eq(sec->afpc[i], prod_b);

      }

    /* Free data */
    cplx_vfree(old_a);
    cplx_vfree(old_b);

    //    for(i = 0; i < s->n; i++)
    //      {
    //        MPS_DEBUG_CPLX(s, sec->afpc[i], "sec->afpc[%d]", i);
    //        MPS_DEBUG_CPLX(s, sec->bfpc[i], "sec->bfpc[%d]", i);
    //      }

    mps_secular_fstart(s, s->n, 0, 0, 0, s->eps_out);

    break;

    /* If this is the DPE phase regenerate DPE coefficients */
  case dpe_phase:

    /* Allocate old_a and old_b */
    old_da = cdpe_valloc(s->n);
    old_db = cdpe_valloc(s->n);

    /* Copy the old coefficients, and set the new
     * b_i with the current roots approximations. */
    for (i = 0; i < s->n; i++)
      {
        cdpe_set(old_da[i], sec->adpc[i]);
        cdpe_set(old_db[i], sec->bdpc[i]);
        cdpe_set(sec->bdpc[i], s->droot[i]);
      }

    /* Compute the new a_i */
    for (i = 0; i < s->n; i++)
      {
        cdpe_t prod_old_b, prod_b, sec_ev;
        cdpe_t ctmp, btmp;
        cdpe_set(prod_old_b, cdpe_one);
        cdpe_set(prod_b, cdpe_one);
        cdpe_set(sec_ev, cdpe_zero);

        for (j = 0; j < sec->n; j++)
          {
            /* Compute 1 / (b - old_b) */
            cdpe_sub(btmp, sec->bdpc[i], old_db[j]);
            cdpe_inv(ctmp, btmp);

            /* Add a_j / (b_i - old_b_j) to sec_ev */
            cdpe_mul_eq(ctmp, old_da[j]);
            cdpe_add_eq(sec_ev, ctmp);

            /* Multiply prod_b for
             * b_i - b_j if i \neq j and prod_old_b
             * for b_i - old_b_i.  */
            cdpe_mul_eq(prod_b, btmp);
            if (i != j)
              {
                cdpe_sub(ctmp, sec->bdpc[i], sec->bdpc[j]);
                cdpe_mul_eq(prod_old_b, ctmp);
              }
          }

        /* Compute the new a_i as sec_ev * prod_old_b / prod_b */
        cdpe_sub_eq(sec_ev, cdpe_one);
        cdpe_mul(sec->adpc[i], sec_ev, prod_old_b);
        cdpe_div_eq(sec->adpc[i], prod_b);

        MPS_DEBUG_CDPE(s, sec_ev, "sec_ev");
      }

    /* Free data */
    cdpe_vfree(old_da);
    cdpe_vfree(old_db);

    /* Debug new coefficients found */
    for (i = 0; i < s->n; i++)
      {
        MPS_DEBUG_CDPE(s, sec->adpc[i], "sec->adpc[%d]", i);
      }

    mps_secular_dstart(s, s->n, 0, (__rdpe_struct *) rdpe_zero,
        (__rdpe_struct *) rdpe_zero, s->eps_out);

    break;

  case mp_phase:
    /* Allocate old_a and old_b */
    old_ma = mpc_valloc(s->n);
    old_mb = mpc_valloc(s->n);

    mpc_vinit2(old_ma, s->n, s->mpwp);
    mpc_vinit2(old_mb, s->n, s->mpwp);

    /* Copy the old coefficients, and set the new
     * b_i with the current roots approximations. */
    for (i = 0; i < s->n; i++)
      {
        mpc_set(old_ma[i], sec->ampc[i]);
        mpc_set(old_mb[i], sec->bmpc[i]);
        mpc_set(sec->bmpc[i], s->mroot[i]);
      }

    /* Compute the new a_i */
    for (i = 0; i < s->n; i++)
      {
        mpc_set_ui(prod_old_b, 1, 0);
        mpc_set_ui(prod_b, 1, 0);
        mpc_set_ui(sec_ev, 0, 0);

        for (j = 0; j < sec->n; j++)
          {
            /* Compute 1 / (b_i - old_b_j) */
            mpc_sub(btmp, sec->bmpc[i], old_mb[j]);

            /* If b - old_b is zero, simplify the computation */
            if (mpc_eq_zero(btmp))
              {
                for (i = 0; i < s->n; i++)
                  {
                    mpc_set(sec->bmpc[i], old_mb[i]);
                  }

                MPS_DEBUG(s, "Cannot regenerate coefficients, reusing old ones.")
                goto regenerate_m_exit;
              }

            mpc_inv(ctmp, btmp);

            /* Add a_j / (b_i - old_b_j) to sec_ev */
            mpc_mul_eq(ctmp, old_ma[j]);
            mpc_add_eq(sec_ev, ctmp);

            /* Multiply prod_b for
             * b_i - b_j if i \neq j and prod_old_b
             * for b_i - old_b_i.  */
            mpc_mul_eq(prod_b, btmp);
            if (i != j)
              {
                mpc_sub(ctmp, sec->bmpc[i], sec->bmpc[j]);
                mpc_mul_eq(prod_old_b, ctmp);
              }
          }

        /* Compute the new a_i as sec_ev * prod_old_b / prod_b */
        mpc_sub_eq_ui(sec_ev, 1, 0);
        mpc_mul(sec->ampc[i], sec_ev, prod_old_b);
        mpc_div_eq(sec->ampc[i], prod_b);

      }

    regenerate_m_exit:

    /* Free data */
    mpc_vclear(old_ma, s->n);
    mpc_vclear(old_mb, s->n);
    mpc_vfree(old_ma);
    mpc_vfree(old_mb);

    for (i = 0; i < s->n; i++)
      {
        MPS_DEBUG_MPC(s, 15, sec->ampc[i], "sec->ampc[%d]", i);
        MPS_DEBUG_MPC(s, 15, sec->bmpc[i], "sec->bmpc[%d]", i);
      }

    //    mps_secular_mstart(s, s->n, 0, (__rdpe_struct *) rdpe_zero,
    //        (__rdpe_struct *) rdpe_zero, s->eps_out);

    break;

  default:
    break;

    } /* End of switch (s->lastphase)*/

  /* Finally set radius according to new computed a_i coefficients,
   * if they are convenient   */
  mps_secular_set_radii(s);
}

void
mps_secular_raise_precision(mps_status* s)
{
  MPS_DEBUG_THIS_CALL

  int i;
  s->mpwp *= 2;
  mps_secular_equation* sec = (mps_secular_equation*) s->user_data;
  for (i = 0; i < s->n; i++)
    {
      mpc_set_prec(sec->ampc[i], s->mpwp);
      mpc_set_prec(sec->bmpc[i], s->mpwp);
      mpc_set_prec(s->mroot[i], s->mpwp);
    }
  rdpe_set_2dl(s->mp_epsilon, 1.0, -s->mpwp);
  MPS_DEBUG(s, "Precision is now at %d bits", s->mpwp);
}

/**
 * @brief Prepare data for the iteration in the new phase specified
 * in the second parameter.
 *
 * Note that for now this function is only able to handle switch
 * from floating point phases (i.e. float_phase or dpe_phase) to
 * multiprecision, and not coming back.
 */
void
mps_secular_switch_phase(mps_status* s, mps_phase phase)
{
  MPS_DEBUG_THIS_CALL

  int i = 0;
  mps_secular_equation* sec = (mps_secular_equation*) s->user_data;
  if (phase == mp_phase)
    {
      s->mpwp = DBL_MANT_DIG;
      mps_secular_raise_precision(s);
      switch (s->lastphase)
        {
      case float_phase:
        /* Copy the approximated roots and the
         * secular equation coefficients */
        for (i = 0; i < s->n; i++)
          {
            mpc_set_cplx(s->mroot[i], s->froot[i]);
            mpc_set_cplx(sec->ampc[i], sec->afpc[i]);
            mpc_set_cplx(sec->bmpc[i], sec->bfpc[i]);
            rdpe_set_d(s->drad[i], s->frad[i]);
          }
        break;

      case dpe_phase:
        /* Copy the coefficients and the approximated
         * roots into the multiprecision values    */
        for (i = 0; i < s->n; i++)
          {
            mpc_set_cdpe(s->mroot[i], s->droot[i]);
            mpc_set_cdpe(sec->ampc[i], sec->adpc[i]);
            mpc_set_cdpe(sec->bmpc[i], sec->bdpc[i]);
          }

      default:
        break;

        }

      /* Set lastphase to mp_phase */
      s->lastphase = mp_phase;

      /* Set epsilon */
      rdpe_set_2dl(s->mp_epsilon, 1.0, -s->mpwp + 1);
    }
  else
    {
      fprintf(stderr, "mps_secular_switch_phase is only able to manage\n"
        "switches from float_phase or dpe_phase to mp_phase. Aborting.");
      exit(EXIT_FAILURE);
    }
}

/*
 * @brief Check if iterations can terminate
 */
mps_boolean
mps_secular_ga_check_stop(mps_status* s)
{
  double frad = pow(10, -s->prec_out);
  rdpe_t drad;
  rdpe_t max_rad;
  int i;
  rdpe_set_2dl(drad, 1.0, -s->prec_out);

  if (frad == 0 && s->lastphase == float_phase)
    return false;

  rdpe_set(max_rad, drad);

  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
        {
      case float_phase:
        if (s->frad[i] > frad)
          return false;
      case mp_phase:
        if (rdpe_gt(s->drad[i], drad))
          {
            if (rdpe_gt(s->drad[i], max_rad))
              rdpe_set(max_rad, s->drad[i]);
            // return false;
          }
        break;
      default:
        break;
        }
    }

  if (rdpe_eq(max_rad, drad))
    return true;

  MPS_DEBUG_RDPE(s, max_rad, "Maximum rad is too big to exit: max_rad")
  return false;
}

/**
 * @brief Update radii of the roots according to the coefficients
 * of the secular equation in this moment, if they are better of
 * the radii present now.
 */
void
mps_secular_set_radii(mps_status* s)
{
  MPS_DEBUG_THIS_CALL

  int i;
  mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

  /* Select right computation based on the phase we are in
   * right now   */
  switch (s->lastphase)
    {

  case float_phase:
    {
      /* Floating point implementation */
      double rad, total_rad = 0;

      /* Compute total radius as \sum_i |sec->afpc[i]| */
      for (i = 0; i < s->n; i++)
        total_rad += cplx_mod(sec->afpc[i]);

      /* Check if the Gerschgorin's radii are more convenient */
      for (i = 0; i < s->n; i++)
        {
          /* TODO: Use the guaranteed computation */
          rad = s->n * cplx_mod(sec->afpc[i]);
          if (rad > total_rad)
            rad = total_rad;

          if (rad < s->frad[i])
            {
              MPS_DEBUG(s, "Setting s->frad[%d] to %e (thanks to Gerschgorin)", i, rad);
              s->frad[i] = rad;
            }
        }
    }
    break;
  case dpe_phase:
  case mp_phase:
    {
      /* DPE and multiprecision implementation */
      rdpe_t rad, total_rad, rtmp;
      cdpe_t ctmp;
      rdpe_set(total_rad, rdpe_zero);

      /* Compute total radius as \sum_i |sec->afpc[i]| */
      for (i = 0; i < s->n; i++)
        {
          if (s->lastphase == mp_phase)
            {
              mpc_get_cdpe(ctmp, sec->ampc[i]);
              cdpe_mod(rtmp, ctmp);
            }
          else
            /* We are in the DPE phase */
            cdpe_mod(rtmp, sec->adpc[i]);

          rdpe_add_eq(total_rad, rtmp);
        }

      /* Check if the Gerschgorin's radii are more convenient */
      for (i = 0; i < s->n; i++)
        {
          /* TODO: Use the guaranteed computation */
          if (s->lastphase == mp_phase)
            {
              mpc_get_cdpe(ctmp, sec->ampc[i]);
              cdpe_mod(rad, ctmp);
            }
          else
            /* We are in the DPE phase */
            cdpe_mod(rad, sec->adpc[i]);

          /* Check which radius is smaller */
          rdpe_mul_eq_d(rad, (double) s->n);
          if (rdpe_gt(rad, total_rad))
            rdpe_set(rad, total_rad);

          /* If the radius is convenient set it */
          //          MPS_DEBUG_RDPE(s, s->drad[i], "Old s->drad[%d]", i)
          //          MPS_DEBUG_RDPE(s, rad, "New Gerschgorin proposal ")
          if (rdpe_lt(rad, s->drad[i]))
            {
              MPS_DEBUG_RDPE(s, rad, "Setting s->drad[%d] to rad", i);
              rdpe_set(s->drad[i], rad);
            }
        }
    }
    break;

  default:
    break;
    }
}

/**
 * @brief MPSolve main function for the secular equation solving
 * using Gemignani's approach.
 */
void
mps_secular_ga_mpsolve(mps_status* s, mps_phase phase)
{
  int roots_computed = 0;
  int iteration_per_packet = 10;
  int i;

  for (i = 0; i < s->n; i++)
    s->frad[i] = DBL_MAX;

  /* Set initial cluster structure as no cluster structure. */
  mps_cluster_reset(s);

  /* Set phase */
  s->lastphase = phase;

  /* Select initial approximations using the custom secular
   * routine and based on the phase selected by the user. */
  switch (s->lastphase)
    {
  case float_phase:
    mps_secular_fstart(s, s->n, 0, 0.0, 0.0, s->eps_out);
    break;

  case dpe_phase:
    mps_secular_dstart(s, s->n, 0, (__rdpe_struct *) rdpe_zero,
        (__rdpe_struct *) rdpe_zero, s->eps_out);
    break;

  case mp_phase:
    mps_secular_mstart(s, s->n, 0, (__rdpe_struct *) rdpe_zero,
        (__rdpe_struct *) rdpe_zero, s->eps_out);
    break;

  default:
    break;
    }

  /* Set initial radius */
  mps_secular_set_radii(s);

  /* Cycle until approximated */
  do
    {
      /* Perform an iteration of floating point Aberth method */
      switch (s->lastphase)
        {
      case float_phase:
        roots_computed = mps_secular_ga_fiterate(s, iteration_per_packet);
        MPS_DEBUG(s, "%d roots were computed", roots_computed)
        break;

      case dpe_phase:
        roots_computed = mps_secular_ga_diterate(s, iteration_per_packet);
        MPS_DEBUG(s, "%d roots were computed", roots_computed)
        break;

      case mp_phase:
        roots_computed = mps_secular_ga_miterate(s, iteration_per_packet);
        MPS_DEBUG(s, "%d roots were computed", roots_computed)
        break;

      default:
        break;
        }

      /* Check if all roots were approximated with the
       * given input precision                      */
      if (mps_secular_ga_check_stop(s))
        return;

      /* Check if it's time to abandon floating point to enter
       * the multiprecision phase */
      if (s->lastphase != mp_phase && (roots_computed == s->n))
        {
          MPS_DEBUG(s, "Switching to multiprecision phase")
          mps_secular_switch_phase(s, mp_phase);

          /* Regenerate coefficients is able to understand the type
           * of data that we are treating, so no switch is necessary
           * in here. */
          mps_secular_ga_regenerate_coefficients(s);
        }
      else if (s->lastphase == mp_phase)
        {
          /* If all the roots were approximated and we are in the multiprecision
           * phase then it's time to increase the precision, or stop if enough
           * precision has been reached. */
          mps_secular_raise_precision(s);

          /* Regenerate coefficients to accelerate convergence. */
          mps_secular_ga_regenerate_coefficients(s);
        }
    }
  while (!mps_secular_ga_check_stop(s));
}
