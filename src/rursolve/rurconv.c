/***********************************************************
**       Multiprecision Polynomial Solver (MPSolve)       **
**              Version 2.1, september 1999               **
**                                                        **
**                      Written by                        **
**       Dario Andrea Bini and Giuseppe Fiorentino        **
**       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
**                                                        **
** (C) 1999, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

/* forward declarations */
void update (int poly, int deg);
void prescan (void);
void scan (void);
void skip (void);
void skipnum (void);
void copynum (void);
void error (void);

/* global variables */
static int polydeg[100];
static int polymon[100];
static int degmax;

FILE *fp = NULL;
int poly;
char c;

/* main program */
int
main (int argc, char *argv[])
{
  char *filename;

  if (argc != 2)
    {
      fprintf (stderr, "Input file missing\n");
      exit (1);
    }
  filename = argv[1];

  if ((fp = fopen (filename, "r")) == NULL)
    {
      fprintf (stderr, "Cannot open input file\n");
      exit (1);
    }
  prescan ();
  rewind (fp);
  scan ();

  fclose (fp);
  return 0;
}

void
prescan ()
{
  int var, prec, deg;

  if (!fscanf (fp, "%d", &var))
    error ();
  if (!fscanf (fp, "%d", &prec))
    error ();

  skip ();
  if (c != '{')
    error ();
  skip ();                      /* eat '{' */
  for (poly = 0; c == '{'; poly++)
    {
      do
        {
          if (!fscanf (fp, "%d", &deg))
            error ();
          update (poly, deg);
          polymon[poly]++;
          skipnum ();
          if (c != ',' && c != '}')
            error ();
        }
      while (c == ',');
      skip ();                  /* eat '}' */
      if (c != ',' && c != '}')
        error ();
      if (c == ',')
        skip ();                /* eat ',' */
    }                           /* for */
  if (c != '}')
    error ();
}

void
scan ()
{
  int var, prec, deg, mon;

  fscanf (fp, "%d", &var);
  if (var != poly - 2)
    error ();
  fscanf (fp, "%d", &prec);
  printf ("%d\n%d\n%d\n", var, degmax, prec);

  skip ();                      /* eat first "{" */
  skip ();                      /* eat '{' */
  for (poly = 0; poly <= var + 1; poly++)
    {
      printf ("\nsri\n");
      if (poly == 0)
        printf ("0\n");
      printf ("%d\n%d\n", polydeg[poly], polymon[poly]);
      for (mon = 0; mon < polymon[poly]; mon++)
        {
          fscanf (fp, "%d", &deg);
          printf (" %d\n", deg);
          copynum ();
        }
      skip ();                  /* eat '}' */
      if (c == ',')
        skip ();                /* eat ',' */
    }
}

void
update (int poly, int deg)
{
  if (polydeg[poly] < deg)
    polydeg[poly] = deg;
  if (deg > degmax && poly > 0)
    degmax = deg;
}

void
skip ()
{
  c = fgetc (fp);
  while (isspace (c))
    c = fgetc (fp);
}

void
skipnum ()
{
  skip ();
  if (c == '-' || c == '+')
    c = fgetc (fp);
  if (!isdigit (c))
    error ();
  while (isdigit (c))
    c = fgetc (fp);
  if (!isspace (c) && c != ',' && c != '}')
    error ();
  if (isspace (c))
    skip ();
}

void
copynum ()
{
  skip ();
  if (c == '-' || c == '+')
    {
      putchar (c);
      c = fgetc (fp);
    }
  while (isdigit (c))
    {
      putchar (c);
      c = fgetc (fp);
    }
  putchar ('\n');
  skip ();
}

void
error ()
{
  fprintf (stderr, "Bad input format\n");
  exit (1);
}
