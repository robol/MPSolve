/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <stdarg.h>
#include <string.h>
#include <mps/mps.h>
#include <ctype.h>
#include <locale.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef __WINDOWS
#include <unistd.h>
#else
#include <io.h>
#endif

#define ISZERO -1

/*********************************************************
*      SUBROUTINE READROOTS                              *
*********************************************************/
void
mps_readroots (mps_context * s)
{
  long digits;
  int i, read_elements;

  if (s->DOLOG)
    fprintf (s->logstr, "Reading roots...\n");

  read_elements = fscanf (s->rtstr, "%ld", &digits);
  if (!read_elements)
    {
      mps_error (s, "Error while reading roots, aborting.");
    }

  /* precision setup code goes here */

  for (i = 0; i < s->n; i++)
    mpc_inp_str_u (s->root[i]->mvalue, s->rtstr, 10);
}

/*********************************************************
*      SUBROUTINE COUNTROOTS                             *
*********************************************************/
void
mps_countroots (mps_context * s)
{
  int k;

  if (s->DOLOG)
    fprintf (s->logstr, "Counting roots...\n");

  s->count[0] = s->count[1] = s->count[2] = 0;

  for (k = 0; k < s->n; k++)
    switch (s->root[k]->inclusion)
      {
      case MPS_ROOT_INCLUSION_IN:
        s->count[0]++;
        break;
      case MPS_ROOT_INCLUSION_OUT:
        s->count[1]++;
        break;
      default:
        s->count[2]++;
        break;
      }

  if (s->output_config->search_set == MPS_SEARCH_SET_UNITARY_DISC_COMPL)
    s->count[1] += s->zero_roots;
  else
    s->count[0] += s->zero_roots;
}

/*********************************************************
*      SUBROUTINE OUTCOUNT                               *
*********************************************************/
void
mps_outcount (mps_context * s)
{
  mps_countroots (s);

  fprintf (s->outstr, "%d roots are inside;\n", s->count[0]);
  fprintf (s->outstr, "%d roots are outside;\n", s->count[1]);
  fprintf (s->outstr, "%d roots are uncertain.\n", s->count[2]);
  if (s->DOLOG)
    {
      fprintf (s->logstr, "%d roots are inside;\n", s->count[0]);
      fprintf (s->logstr, "%d roots are outside;\n", s->count[1]);
      fprintf (s->logstr, "%d roots are uncertain.\n", s->count[2]);
    }
}

/*********************************************************
*      SUBROUTINE OUTFLOAT                               *
*********************************************************/
void
mps_outfloat (mps_context * s, mpf_t f, rdpe_t rad, long out_digit,
              mps_boolean sign)
{
  mpf_t t;
  rdpe_t r, ro;
  double d;
  long l, digit, true_digit;

  if (s->output_config->format == MPS_OUTPUT_FORMAT_FULL)
    {
      mpf_init2 (t, mpf_get_prec (f));
      mpf_set (t, f);
      mpf_out_str (s->outstr, 10, 0, t);
      mpf_clear (t);
      return;
    }

  mpf_init2 (t, s->output_config->prec);

  mpf_get_rdpe (ro, f);
  if (s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT ||
      s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    rdpe_out_str_u (s->outstr, ro);
  else
    {
      rdpe_abs_eq (ro);
      if (rdpe_ne (ro, rdpe_zero))
        rdpe_div (r, rad, ro);
      else
        rdpe_set_d (r, 1.0e-10);
      digit = (long) (-rdpe_log10 (r) - 0.5);
      if (digit <= 0)
        {
          rdpe_get_dl (&d, &l, ro);
          fprintf (s->outstr, "0.e%ld", l);
        }
      else
        {
          true_digit = (long) (LOG10_2 * mpf_get_prec (f));
          true_digit = MIN (digit, true_digit);
          true_digit = MIN (true_digit, out_digit);
          if (sign)
            mpf_set (t, f);
          else
            mpf_abs (t, f);
          mpf_out_str (s->outstr, 10, true_digit, t);
        }
    }

  mpf_clear (t);
}

/*********************************************************
*      SUBROUTINE OUTROOT                                *
*********************************************************/
void
mps_outroot (mps_context * s, int i, int num)
{
  long out_digit;

  out_digit = (long) (LOG10_2 * s->output_config->prec) + 10;

  /* print format header */
  switch (s->output_config->format)
    {
    case MPS_OUTPUT_FORMAT_COMPACT:
    case MPS_OUTPUT_FORMAT_FULL:
      fprintf (s->outstr, "(");
      break;
    case MPS_OUTPUT_FORMAT_VERBOSE:
      fprintf (s->outstr, "Root(%d) = ", num);
      break;
    default:
      break;
    }

  /* print real part */
  if (i == ISZERO || s->root[i]->attrs == MPS_ROOT_ATTRS_IMAG)
    fprintf (s->outstr, "0");
  else
    mps_outfloat (s, mpc_Re (s->root[i]->mvalue), s->root[i]->drad, out_digit, true);

  /* print format middle part */
  switch (s->output_config->format)
    {
    case MPS_OUTPUT_FORMAT_BARE:
      fprintf (s->outstr, " ");
      break;
    case MPS_OUTPUT_FORMAT_GNUPLOT:
    case MPS_OUTPUT_FORMAT_GNUPLOT_FULL:
      fprintf (s->outstr, "\t");
      break;
    case MPS_OUTPUT_FORMAT_COMPACT:
    case MPS_OUTPUT_FORMAT_FULL:
      fprintf (s->outstr, ", ");
      break;
    case MPS_OUTPUT_FORMAT_VERBOSE:
      if (i == ISZERO || mpf_sgn (mpc_Im (s->root[i]->mvalue)) >= 0)
        fprintf (s->outstr, " + I * ");
      else
        fprintf (s->outstr, " - I * ");
      break;
    default:
      break;
    }

  /* print imaginary part */
  if (i == ISZERO || s->root[i]->attrs == MPS_ROOT_ATTRS_REAL)
    fprintf (s->outstr, "0");
  else
    mps_outfloat (s, mpc_Im (s->root[i]->mvalue), s->root[i]->drad, out_digit,
                  s->output_config->format != MPS_OUTPUT_FORMAT_VERBOSE);

  /* If the output format is GNUPLOT_FORMAT_FULL, print out also the radius */
  if (s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      fprintf (s->outstr, "\t");
      rdpe_out_str_u (s->outstr, s->root[i]->drad);
      fprintf (s->outstr, "\t");
      rdpe_out_str_u (s->outstr, s->root[i]->drad);
    }

  /* print format ending */
  switch (s->output_config->format)
    {
    case MPS_OUTPUT_FORMAT_COMPACT:
      fprintf (s->outstr, ")");
      break;
    case MPS_OUTPUT_FORMAT_FULL:
      fprintf (s->outstr, ")\n");
      if (i != ISZERO)
        {
          rdpe_outln_str (s->outstr, s->root[i]->drad);
          fprintf (s->outstr, "Status: %s, %s, %s\n", 
                   MPS_ROOT_STATUS_TO_STRING (s->root[i]->status),
                   MPS_ROOT_ATTRS_TO_STRING (s->root[i]->attrs),
                   MPS_ROOT_INCLUSION_TO_STRING (s->root[i]->inclusion));
        }
      else
        fprintf (s->outstr, " 0\n ---\n");
      break;
    default:
      break;
    }
  fprintf (s->outstr, "\n");

  /* debug info */
  if (s->DOLOG)
    {
      if (i == ISZERO)
        fprintf (s->logstr, "zero root %-4d = 0", num);
      else
        {
          fprintf (s->logstr, "Root %-4d = ", i);
          mpc_out_str_2 (s->logstr, 10, 0, 0, s->root[i]->mvalue);
          fprintf (s->logstr, "\n");
          fprintf (s->logstr, "  Radius = ");
          rdpe_outln_str (s->logstr, s->root[i]->drad);
          fprintf (s->logstr, "  Prec = %ld\n",
                   (long) (mpc_get_prec (s->root[i]->mvalue) / LOG2_10));
          fprintf (s->logstr, "  Approximation = %s\n", 
                   MPS_ROOT_STATUS_TO_STRING (s->root[i]->status));
          fprintf (s->logstr, "  Attributes = %s\n",
                   MPS_ROOT_ATTRS_TO_STRING (s->root[i]->attrs));
          fprintf (s->logstr, "  Inclusion = %s\n",
                   MPS_ROOT_INCLUSION_TO_STRING (s->root[i]->inclusion));
          fprintf (s->logstr, "--------------------\n");
        }
    }
}

/*********************************************************
*      SUBROUTINE OUTPUT                                 *
*********************************************************/
void
mps_output (mps_context * s)
{
  int i, ind, num = 0;

  if (s->DOLOG)
    fprintf (s->logstr, "--------------------\n");

  if (s->output_config->format != MPS_OUTPUT_FORMAT_GNUPLOT && 
      s->output_config->format != MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      if (s->over_max)
        {
          mps_warn (s, "Warning: Input precision has been reached during computation, "
                    "so not all the required digits may have been computed.");
        }
    }

  /* Start with plotting instructions in the case of 
   * MPS_OUTPUT_GNUPLOT_FULL, so the output can be
   * piped directly to gnuplot */
  if (s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      fprintf (s->outstr, "# MPSolve output for GNUPLOT\n");
      fprintf (s->outstr, "# Make user that this output is piped into gnuplot using a command like\n");
      fprintf (s->outstr, "# mpsolve -Ogf | gnuplot \n");
      fprintf (s->outstr, "set pointsize 0.3\n");
      fprintf (s->outstr, "plot '-' title 'Computed roots' with %s\n", s->gnuplot_format);
    }

  if (s->output_config->goal == MPS_OUTPUT_GOAL_COUNT)
    mps_outcount (s);
  else
    {
      if (s->output_config->search_set != MPS_SEARCH_SET_UNITARY_DISC_COMPL)
        for (i = 0; i < s->zero_roots; i++)
            mps_outroot (s, ISZERO, num++);
      for (ind = 0; ind < s->n; ind++)
        {
          i = s->order[ind];
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_OUT)
            continue;
          mps_outroot (s, i, num++);
        }
    }

  if (s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      fprintf (s->outstr, "e\n");
      fprintf (s->outstr, "pause mouse close\n");
      fprintf (s->outstr, "# End of MPSolve GNUPLOT output. If you are seeing this maybe\n");
      fprintf (s->outstr, "# you forgot to pipe the ***solve command into gnuplot?\n");
    }
}

/*********************************************************
*      SUBROUTINE COPY_ROOTS                             *
*********************************************************/
void
mps_copy_roots (mps_context * s)
{
  int i;

  MPS_DEBUG_THIS_CALL;

  switch (s->lastphase)
    {
    case no_phase:
      mps_error (s, "Nothing to copy");
      break;

    case float_phase:
      if (s->DOSORT)
        mps_fsort (s);
      for (i = 0; i < s->n; i++)
        {
          mpc_set_prec (s->root[i]->mvalue, DBL_MANT_DIG);
          mpc_set_cplx (s->root[i]->mvalue, s->root[i]->fvalue);
          rdpe_set_d (s->root[i]->drad, s->root[i]->frad);
        }
      break;

    case dpe_phase:
      if (s->DOSORT)
        mps_dsort (s);
      for (i = 0; i < s->n; i++)
        {
          mpc_set_prec (s->root[i]->mvalue, DBL_MANT_DIG);
          mpc_set_cdpe (s->root[i]->mvalue, s->root[i]->dvalue);
        }
      break;

    case mp_phase:
      if (s->DOSORT)
        mps_msort (s);
      break;

    }
}

/*************************************************************
 *                     SUBROUTINE DUMP                       *
 *************************************************************/
void
mps_dump (mps_context * s)
{
  int i;
  FILE * dmpstr = s->logstr;

  MPS_DEBUG (s, "Dumping the approximations:");

  /* output current status */
  /* fprintf (dmpstr, */
  /*          "Phase=%d, In=%d, Out=%d, Uncertain=%d, Zero=%d, Clusters=%ld\n", */
  /*          s->lastphase, s->count[0], s->count[1], s->count[2], s->zero_roots, */
  /*          s->clusterization->n); */

  MPS_DEBUG (s, "Phase = %s, In = %d, Out = %d, Uncertain = %d, Zero = %d, Cluster = %ld",
             MPS_PHASE_TO_STRING (s->lastphase), s->count[0], s->count[1], s->count[2],
             s->zero_roots, s->clusterization->n);

  /* output current approximations */
  /* fprintf (dmpstr, "\nCurrent approximations:\n"); */
  MPS_DEBUG (s, "Current approximations:");
  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
        {
        case no_phase:
        case float_phase:
          MPS_DEBUG_CPLX (s, s->root[i]->fvalue, "Approximation  %4d", i);
          break;

        case dpe_phase:
          MPS_DEBUG_CDPE (s, s->root[i]->dvalue, "Approximation  %4d", i);
          break;

        case mp_phase:
          MPS_DEBUG_MPC (s, s->mpwp, s->root[i]->mvalue, "Approximation  %4d", i);
          break;
        }
    }

  /* output radii */
  MPS_DEBUG (s, "Current radii:");
  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
        {
        case no_phase:
        case float_phase:
          MPS_DEBUG (s, "Radius of root %4d = %e", i, s->root[i]->frad);
          break;

        case dpe_phase:
        case mp_phase:
          MPS_DEBUG_RDPE (s, s->root[i]->drad, "Radius of root %4d", i);
          break;
        }
    }

  MPS_DEBUG (s, " ");
  mps_dump_status (s, dmpstr);
}

/**
 * @brief Dump cluster structure to <code>outstr</code>.
 *
 * @param s the mps_context struct pointer.
 * @param outstr The output stream where the cluster structure
 *  will be dumped.
 */
void
mps_dump_cluster_structure (mps_context * s, FILE * outstr)
{
  fprintf (outstr,
           "    MPS_DUMP_CLUSTER_STRUCTURE: Dumping cluster structure\n");

  mps_cluster_item * cluster_item;
  mps_cluster * cluster;
  mps_root * root;

  for (cluster_item = s->clusterization->first; cluster_item != NULL;
       cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;
      fprintf (outstr, "     Cluster contains %ld roots:\n", cluster->n);

      /* Dump cluster roots, but not more than 15 for line, to make
       * the output readable. */
      int j = 0;
      for (root = cluster->first; root != NULL; root = root->next)
        {
          /* Go to a newlint if 15 roots are printed out */
          if (j % 15 == 0)
            {
              fprintf (outstr, "\n       ");
            }

          fprintf (outstr, " %4ld", root->k);
          j++;
        }

      /* Make space untile the next cluster */
      fprintf (outstr, "\n\n");
    }
}

/**
 * @brief Dump status of all the root approximations
 */
void
mps_dump_status (mps_context * s, FILE * outstr)
{
  int i;
  MPS_DEBUG (s, "              Approximation              Attributes       Inclusion");
  for (i = 0; i < s->n; i++)
    {
      MPS_DEBUG (s, "Status  %4d: %-25s  %-15s  %-15s", i,
                 MPS_ROOT_STATUS_TO_STRING (s->root[i]->status),
                 MPS_ROOT_ATTRS_TO_STRING  (s->root[i]->attrs), 
                 MPS_ROOT_INCLUSION_TO_STRING (s->root[i]->inclusion));
    }
}

/*************************************************************
 *                     SUBROUTINE WARN                       *
 *************************************************************/
void
mps_warn (mps_context * st, char *s)
{
  if (st->DOWARN)
    {
      if (s[strlen (s)] == '\n')
        {
          fprintf (stderr, "%s", s);
        }
      else
        {
          fprintf (stderr, "%s\n", s);
        }
    }
}

/**
 * @brief Check if the file descriptor associated to stream
 * is bounded to a tty.
 *
 * @param stream the stream to check
 */
mps_boolean
mps_is_a_tty (FILE * stream)
{
#ifndef __WINDOWS
  return isatty (fileno(stream));
#else
  return _isatty (_fileno (stream));
#endif
}

/* Prepare a implicit definition if not provided by the compiler but
 * available as a non-conformant extension. */
#ifndef vsnprintf
#ifdef HAVE_VSNPRINTF
int snprintf(char *str, size_t size, const char *format, ...);
#endif
#endif

/*************************************************************
 *                     SUBROUTINE VAERROR                    *
 *************************************************************/
void
mps_error (mps_context * s, const char * format, ...)
{
  va_list ap;
  int buffer_size = 32;
  int missing_characters = 0;

  va_start (ap, format);

  s->error_state = true;
  if (s->last_error == NULL)
    s->last_error = mps_newv (char, buffer_size);

  /* Measure space needed for the string, if our initial guess for the space neede
   * is not enough */
  while ((missing_characters = vsnprintf (s->last_error, buffer_size, format, ap)) > buffer_size)
  {
    buffer_size += missing_characters + 1;
    s->last_error = mps_realloc (s->last_error, buffer_size);
  }

  va_end (ap);
}

void
mps_print_errors (mps_context * s)
{
  if (s->logstr == NULL)
    s->logstr = stderr;

  if (mps_is_a_tty (s->logstr))
    mps_warn (s, "\033[31;1m!\033[0m MPSolve encountered an error:");  /* output error message */
  else
    mps_warn (s, "! MPSolve encountered an error:");

  mps_warn (s, s->last_error);

  /* Dump approximations, but only if they are present */
  if (s->root && s->lastphase)
    mps_dump (s);
}
