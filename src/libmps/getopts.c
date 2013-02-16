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


#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <float.h>
#include <mps/mps.h>

/**
 * @brief Parse command line options in a similar way of getopts.
 *
 * @param opt_ptr A pointer to a previous used mps_opts that can
 * be re-used, or NULL. If this is NULL it will be allocated.
 * When this function return false the opt_ptr is automatically
 * freed.
 * @param argc_ptr The address in memory of the argc variable
 * obtained by the operating system.
 * @param argv_ptr The address in memory of the argv variable
 * obtained by the operating system.
 * @param opt_format The option that <code>mps_getopts()</code>
 * has to look for. It is passed as a string like
 * <code>"ab:n"</code> where the <code>":"</code> means that the preceding
 * character needs an argument, while single characters not followed
 * by <code>":"</code> don't expect one. In this case a correct call
 * to the program would be
 * @code
 *   ./myprogram -a -n -b 7
 * @endcode
 * or even
 * @code
 *   ./myprogam -a
 * @endcode
 *
 * The option are parsed from argc_ptr and argv_ptr and these
 * are modified taking options away (as the program has been
 * called wihtout these options) so after option parsing the program
 * can still parse other options such as filenames or similar, without
 * worrying about options switch.
 */
mps_boolean
mps_getopts (mps_opt ** opt_ptr, int *argc_ptr, char ***argv_ptr,
             const char *opt_format)
{
  char **argv = *argv_ptr;
  int argc = *argc_ptr;
  char *tmp;
  mps_opt *opt;
  int i, l = strlen (opt_format), steps;

  if ((*opt_ptr) == NULL)
    {
      (*opt_ptr) = (mps_opt *) mps_malloc (sizeof (mps_opt));
    }

  opt = *opt_ptr;

  /* Check if there are other arguments to parse */
  if (!argc)
    {
      free (opt);
      return false;
    }

  if (argc == 1)
    {
      free (opt);
      return false;
    }

  /* Scan for right offset */
  steps = 0;
  while (argv[1][0] != '-')
    {
      /* If we get here then argv[1] is a parameter string
       * that must be preserved, so we should put it on the end
       * of argv. */
      tmp = argv[1];
      for (i = 1; i < argc - 1; i++)
        {
          (*argv_ptr)[i] = argv[i + 1];
        }
      (*argv_ptr)[argc - 1] = tmp;
      steps++;

      /* Check if argc permutation was performed */
      if (steps == argc - 1)
        {
          free (opt);
          return false;
        }
    }


  /* Search argv[0][1] in opt_format */
  for (i = 0; i < l; i++)
    {
      if (opt_format[i] == ':')
        continue;
      else if (argv[1][1] == opt_format[i])
        {
          opt->optchar = opt_format[i];

          /* If there is no argument we are done */
          if (i == l || opt_format[i + 1] != ':')
            {
              opt->optvalue = NULL;
              (*argc_ptr)--;
              (*argv_ptr)[1] = argv[0];
              (*argv_ptr)++;
              return true;
            }

          /* Check if the parameter has an optional argument. If that's the
           * case, don't bother if there are no arguments.
           */
          if (i <= l + 2
              && (opt_format[i + 1] == ':' && opt_format[i + 2] == ':'))
            {
              if (argv[1][2] == '\0')
                {
                  opt->optvalue = NULL;
                  (*argc_ptr)--;
                  (*argv_ptr)[1] = argv[0];
                  (*argv_ptr)++;
                  return true;
                }
            }

          /* If the string is not terminated than we should
           * expect to find the parameter attached to it */
          if (argv[1][2] != '\0')
            {
              if (argv[1][2] == '=')
                {
                  opt->optvalue = argv[1] + 3;
                  (*argc_ptr)--;
                  (*argv_ptr)[1] = argv[0];
                  (*argv_ptr)++;
                  return true;
                }
              else
                {
                  opt->optvalue = argv[1] + 2;
                  (*argc_ptr)--;
                  (*argv_ptr)[1] = argv[0];
                  (*argv_ptr)++;
                  return true;
                }
            }

          /* If the parameter should be in argv[1] but is not
           * there return an error */
          if (argc == 2)
            {
              opt->optvalue = NULL;
              (*argc_ptr)--;
              (*argv_ptr)[1] = argv[0];
              (*argv_ptr)++;
              return true;
            }

          /* Otherwise we can set the argument in the struct */
          opt->optvalue = argv[2];
          (*argc_ptr) -= 2;
          (*argv_ptr)[2] = argv[0];
          (*argv_ptr) += 2;
          return true;
        }
    }

  /* Fallback case, we haven't recognized the input */
  opt->optchar = '?';
  (*argc_ptr)--;
  (*argv_ptr)[1] = argv[0];
  (*argv_ptr)++;

  return true;
}
