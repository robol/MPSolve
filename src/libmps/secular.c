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
#include <ctype.h>
#include <mps/mpc.h>

/**
 * @brief Deflate a secular equation lowering the degree of the
 * polynomial that represent it, if that is possible.
 */
void
mps_secular_deflate(mps_status* s, mps_secular_equation* sec)
{
  int i, j, k;
  for (i = 0; i < sec->n; i++)
    {
      for (j = i + 1; j < sec->n; j++)
        {
          if (cdpe_eq(sec->bdpc[i], sec->bdpc[j]))
            {
              cdpe_add_eq(sec->adpc[i], sec->adpc[j]);

              /* Copy other coefficients back of one position */
              for (k = j; k < sec->n - 1; k++)
                {
                  cdpe_set(sec->adpc[j], sec->adpc[j + 1]);
                  cdpe_set(sec->bdpc[j], sec->bdpc[j + 1]);
                }

              /* Decrement number of coefficients */
              sec->n--;
            }
        }
    }
}

/**
 * @brief Raw version of mps_secular_equation_new that only
 * allocate space for the coefficients but relies on the user
 * to fill their values.
 */
mps_secular_equation*
mps_secular_equation_new_raw(mps_status* s, unsigned long int n)
{
  mps_secular_equation* sec = (mps_secular_equation*) malloc(
      sizeof(mps_secular_equation));

  /* Allocate floating point coefficients */
  sec->afpc = cplx_valloc(n);
  sec->bfpc = cplx_valloc(n);
  sec->old_afpc = cplx_valloc(n);
  sec->old_bfpc = cplx_valloc(n);

  /* Allocate complex dpe coefficients of the secular equation */
  sec->adpc = cdpe_valloc(n);
  sec->bdpc = cdpe_valloc(n);
  sec->old_adpc = cdpe_valloc(n);
  sec->old_bdpc = cdpe_valloc(n);

  /* Allocate multiprecision complex coefficients of the secular equation */
  sec->ampc = mpc_valloc(n);
  sec->bmpc = mpc_valloc(n);
  sec->old_ampc = mpc_valloc(n);
  sec->old_bmpc = mpc_valloc(n);



  /* Init multiprecision arrays */
  mpc_vinit(sec->ampc, n);
  mpc_vinit(sec->bmpc, n);
  mpc_vinit(sec->old_ampc, n);
  mpc_vinit(sec->old_bmpc, n);

  sec->n = n;
  return sec;
}

/**
 * @brief Create a new secular equation struct
 */
mps_secular_equation*
mps_secular_equation_new(mps_status* s, cplx_t* afpc, cplx_t* bfpc, unsigned long int n)
{

  int i;

  /* Allocate the space for the new struct */
  mps_secular_equation* sec = mps_secular_equation_new_raw(s, n);

  /* Copy the complex coefficients passed as argument */
  for (i = 0; i < n; i++)
    {
      /* a_i coefficients */
      cplx_set(sec->afpc[i], afpc[i]);

      /* b_i coefficients */
      cplx_set(sec->bfpc[i], bfpc[i]);
    }

  sec->n = n;
  mps_secular_deflate(s, sec);

  for (i = 0; i < sec->n; i++)
    {
      cdpe_init(sec->adpc[i]);
      cdpe_set_x(sec->adpc[i], sec->afpc[i]);

      mpc_set_cplx(sec->ampc[i], sec->afpc[i]);

      cdpe_init(sec->bdpc[i]);
      cdpe_set_x(sec->bdpc[i], sec->bfpc[i]);

      mpc_set_cplx(sec->bmpc[i], sec->bfpc[i]);
    }

  return sec;
}

void
mps_secular_equation_free(mps_secular_equation* s)
{
  /* Free internal data */
  cplx_vfree(s->afpc);
  cplx_vfree(s->bfpc);

  cdpe_vfree(s->adpc);
  cdpe_vfree(s->bdpc);

  mpc_vclear(s->ampc, s->n);
  mpc_vclear(s->bmpc, s->n);

  mpc_vfree(s->ampc);
  mpc_vfree(s->bmpc);

  /* And old coefficients */
  cplx_vfree(s->old_afpc);
  cplx_vfree(s->old_bfpc);
  cdpe_vfree(s->old_adpc);
  cdpe_vfree(s->old_bdpc);
  mpc_vclear(s->old_ampc, s->n);
  mpc_vclear(s->old_bmpc, s->n);
  mpc_vfree(s->old_ampc);
  mpc_vfree(s->old_bmpc);

  /* ...and then release it */
  free(s);
}

void
mps_skip_comments (FILE* input_stream)
{
    char buf;
    while ((buf = fgetc(input_stream)) == '!' || isspace(buf))
        if (buf == '!')
            /* Skip until newline */
            while(fgetc(input_stream) != '\n');
    ungetc(buf, input_stream);
}

/**
 * @brief Parse a line of the input stream that contains the character
 * ';', so should be considered an option line.
 *
 * Valid options, recognized at the moment being are:
 */
mps_flag
mps_parse_option_line (mps_status* s, char* line, size_t length)
{
  char* first_comment;
  char *option;
  char *c_ptr;
  char buf[255];
  size_t real_length;

  if (length > 255)
      mps_error(s, 1, "Maximum line length exceeded (length > 255 while parsing)");

  /* Check if there are comments in this line */
  if ((first_comment = strchr(line, '!')) != NULL)
      real_length = (first_comment - line) / sizeof(char);
  else
      real_length = length;

  /* Get the characters before ';' and strip spaces */
  c_ptr = line;
  while (isspace(*c_ptr) && ((c_ptr < first_comment) || first_comment == NULL))
  {
      c_ptr++;
      real_length--;
  }
  option = c_ptr;
  c_ptr = strchr(option, ';') - 1;
  while (isspace(*--c_ptr) && real_length--);

  /* Now we have the option that is pointed by option and is
   * real_lenght characters long */
  *(c_ptr + 2) = '\0';
  MPS_DEBUG(s, "Parsed option: %s", option);

  /* Detect option about density-sparseness */
  if (strcasecmp(option, "dense") == 0)
      return MPS_FLAG_DENSE;
  if (strcasecmp(option, "sparse") == 0)
      return MPS_FLAG_SPARSE;

  /* Options on types */
  if (strcasecmp(option, "integer") == 0)
      return MPS_FLAG_INTEGER;
  if (strcasecmp(option, "real") == 0)
      return MPS_FLAG_REAL;
  if (strcasecmp(option, "rational") == 0)
      return MPS_FLAG_RATIONAL;

  /* Options on the input type */
  if (strcasecmp(option, "secular") == 0)
      return MPS_FLAG_SECULAR;
  if (strcasecmp(option, "polynomial") == 0)
      return MPS_FLAG_POLYNOMIAL;

  /* ...if not identified */
  return MPS_FLAG_UNDEFINED;
}

mps_secular_equation*
mps_secular_equation_read_from_stream(mps_status* s, FILE* input_stream)
{
  mps_secular_equation* sec;
  int i, r, n;
  mps_boolean parsing_options = true;
  mps_flag flag;
  size_t length;
  mps_input_buffer *buffer;

  /* Create the input buffer */
  buffer = mps_input_buffer_new (input_stream);

  /* Read options, if present */
  mps_skip_comments(input_stream);
  while (parsing_options)
  {
    // r = getline(&line, &length, input_stream);
    mps_input_buffer_readline(buffer);
    if (strchr(buffer->line, ';') == NULL ||
            mps_input_buffer_eof(buffer))
    {
        parsing_options = false;
    }
    else
    {
        flag = mps_parse_option_line(s, buffer->line, length);
    }
  }

  /* Read the number of the coefficients */
  mps_skip_comments(input_stream);
  while (!r || n == 0)
  {
    r = sscanf(buffer->line, "%d", &n);

    if (!r || n == 0)
        mps_input_buffer_readline(buffer);
  }

  if (!r)
      mps_error(s, 1, "Error reading input coefficients of the secular equation.\n");

  MPS_DEBUG(s, "Degree: %d", n)

  /* Read directly the secular equation in DPE, so we don't need
   * to have a fallback case if the coefficients are bigger than
   * what is supported by the standard floating point arithmetic */
  sec = mps_secular_equation_new_raw(s, n);

  for(i = 0; i < n; i++)
    {
      mps_skip_comments(input_stream);
      rdpe_inp_str_flex(cdpe_Re(sec->adpc[i]), input_stream);

      mps_skip_comments(input_stream);
      rdpe_inp_str_flex(cdpe_Im(sec->adpc[i]), input_stream);

      mps_skip_comments(input_stream);
      rdpe_inp_str_flex(cdpe_Re(sec->bdpc[i]), input_stream);

      mps_skip_comments(input_stream);
      rdpe_inp_str_flex(cdpe_Im(sec->bdpc[i]), input_stream);
    }

  /* Deflate input, if identical b_i coefficients are found */
  mps_secular_deflate(s, sec);

  /* Copy coefficients back in other places */
  for(i = 0; i < sec->n; i++)
    {
      mpc_set_cdpe(sec->ampc[i], sec->adpc[i]);
      mpc_set_cdpe(sec->bmpc[i], sec->bdpc[i]);

      /* Get floating points coefficients */
      cdpe_get_x(sec->afpc[i], sec->adpc[i]);
      cdpe_get_x(sec->bfpc[i], sec->bdpc[i]);
    }

  return sec;
}

/**
 * @brief Utility function that save a snapshot of the coefficient
 * of the secular equation in the fields old_* in the
 * <code>mps_secular_equation</code> struct.
 *
 * Those coefficients will then be used when recomputing the new
 * coefficients in <code>mps_secular_ga_regenerate_coefficients()</code>.
 */
void
mps_secular_save_coefficients (mps_status* s, mps_secular_equation* sec)
{
    int i;
    for(i = 0; i < sec->n; i++)
      {
        cplx_set(sec->old_afpc[i], sec->afpc[i]);
        cplx_set(sec->old_bfpc[i], sec->bfpc[i]);

        cdpe_set(sec->old_adpc[i], sec->adpc[i]);
        cdpe_set(sec->old_bdpc[i], sec->bdpc[i]);

        mpc_set(sec->old_ampc[i], sec->ampc[i]);
        mpc_set(sec->old_bmpc[i], sec->bmpc[i]);
      }
}


/**
 * @brief Evaluate secular equation in the point x.
 */
void
mps_secular_evaluate(mps_status* s, cplx_t x, cplx_t sec_ev)
{
  cplx_t ctmp;
  int i;
  mps_secular_equation* sec = (mps_secular_equation*) s->secular_equation;
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
  mps_secular_equation* sec = (mps_secular_equation*) s->secular_equation;
  *which_case = (sec->starting_case == float_phase) ? 'f' : 'd';
}

void
mps_secular_raise_coefficient_precision(mps_status* s, int wp)
{
  MPS_DEBUG_THIS_CALL

  int i;
  mps_secular_equation* sec = (mps_secular_equation*) s->secular_equation;

  if (s->lastphase != mp_phase)
  {
      MPS_DEBUG(s, "We are still in a floting point phase, so copying the coefficients")
      for (i = 0; i < s->n; i++)
      {
          switch (s->lastphase)
          {
            case dpe_phase:
              mpc_set_cdpe(sec->ampc[i], sec->adpc[i]);
              mpc_set_cdpe(sec->bmpc[i], sec->bdpc[i]);
              mpc_set_cdpe(s->mroot[i], s->droot[i]);
              break;
          case float_phase:
              mpc_set_cplx(sec->ampc[i], sec->afpc[i]);
              mpc_set_cplx(sec->bmpc[i], sec->bfpc[i]);
              mpc_set_cplx(s->mroot[i], s->froot[i]);
              break;
          }
      }
  }

  for (i = 0; i < s->n; i++)
    {
      mpc_set_prec(sec->ampc[i], s->mpwp);
      mpc_set_prec(sec->bmpc[i], s->mpwp);
      mpc_set_prec(s->mroot[i], s->mpwp);

      mpc_set_prec(sec->old_ampc[i], s->mpwp);
      mpc_set_prec(sec->old_bmpc[i], s->mpwp);
    }
  rdpe_set_2dl(s->mp_epsilon, 1.0, -s->mpwp);
  MPS_DEBUG(s, "Precision of the coefficients is now at %d bits", s->mpwp);

  if (s->lastphase != mp_phase)
  {
      MPS_DEBUG(s, "We are still in a floting point phase, so copying the coefficients back")
      for (i = 0; i < s->n; i++)
      {
          switch (s->lastphase)
          {
            case dpe_phase:
              mpc_get_cdpe(sec->adpc[i], sec->ampc[i]);
              mpc_get_cdpe(sec->bdpc[i], sec->bmpc[i]);
              break;
          case float_phase:
              mpc_get_cplx(sec->afpc[i], sec->ampc[i]);
              mpc_get_cplx(sec->bfpc[i], sec->bmpc[i]);
              break;
          }
      }
  }
}

void
mps_secular_raise_precision(mps_status* s, int wp)
{
    mps_secular_raise_coefficient_precision(s, wp);
    s->mpwp = wp;
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
  mps_secular_equation* sec = (mps_secular_equation*) s->secular_equation;
  if (phase == mp_phase)
    {
      s->mpwp = DBL_MANT_DIG;
      mps_secular_raise_precision(s, 2 * s->mpwp);
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
  mps_secular_equation* sec = (mps_secular_equation*) s->secular_equation;

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

      /* Compute guaranteed total rad */
      rdpe_mul_d(rtmp, s->mp_epsilon, s->n);
      rdpe_add_eq(total_rad, rtmp);

      /* Check if the Gerschgorin's radii are more convenient */
      for (i = 0; i < s->n; i++)
        {
          /* TODO: Use the guaranteed computation */
          if (s->lastphase == mp_phase)
            {
              mpc_get_cdpe(ctmp, sec->ampc[i]);
              cdpe_mod(rad, ctmp);
              rdpe_add_eq(rad, s->mp_epsilon);
            }
          else
            {
              /* We are in the DPE phase */
              cdpe_mod(rad, sec->adpc[i]);
              rdpe_add_eq(rad, s->mp_epsilon);
            }

          /* Check which radius is smaller (here guaranteed radius is
           * computed). */
          rdpe_mul_d(rtmp, s->mp_epsilon, 9 * s->n);
          rdpe_add_eq(rtmp, rdpe_one);
          rdpe_mul_eq_d(rad, (double) s->n);
          rdpe_mul_eq(rad, rtmp);

          if (rdpe_gt(rad, total_rad))
            rdpe_set(rad, total_rad);

          /* If the radius is convenient set it */
          if (rdpe_lt(rad, s->drad[i]))
            {
              rdpe_set(s->drad[i], rad);
            }
        }
    }
    break;

  default:
    break;
    }
}
