/***********************************************************
**        Memory Manager for mpf and mpc Temporaries      **
**                      Version 1.0                       **
**                                                        **
**             Written by Giuseppe Fiorentino             **
**                 (fiorent@dm.unipi.it)                  **
**                                                        **
***********************************************************/

#include <mps/mps_mptemp.h>       /* NOMPTEMP may be defined here */
#ifndef NOMPTEMP

#include <stdlib.h>

#ifndef MPTEMP_DEBUG
#define MPTEMP_DEBUG 1
#endif

#ifndef MPTEMP_CHECK
#define MPTEMP_CHECK 1
#endif

#ifndef MPTEMP_RENEW
#define MPTEMP_RENEW 1
#endif

/* number of temporary variables */
#define MINBUFFSIZE 50

/***********************************************************
**              mpf_t temporaries                         **
***********************************************************/

static mpf_t *mpfbuff = NULL;
static long  *mpfprec = NULL;
static int mpfsize = 0;
static int mpfstkp = 0;
static int mpfallp = 0;

tmpf_t
gettmpf(long prec)
{
  if (mpfstkp == mpfallp) {
    if (mpfallp == mpfsize) {
      if (mpfsize) {
	mpfbuff = (mpf_t *) realloc(mpfbuff, 2 * mpfsize * sizeof(mpf_t));
	mpfprec = (long *) realloc(mpfprec, 2 * mpfsize * sizeof(long));
	mpfsize *= 2;
      } else {
	mpfbuff = (mpf_t *) malloc(MINBUFFSIZE * sizeof(mpf_t));
	mpfprec = (long *) malloc(MINBUFFSIZE * sizeof(long));
	mpfsize = MINBUFFSIZE;
      }
    }
    mpf_init2(mpfbuff[mpfallp], prec);
    prec = mpf_get_prec(mpfbuff[mpfallp]);
    mpfprec[mpfallp] = prec;
    mpfallp++;
  }

  if (mpfprec[mpfstkp] >= prec) 
    mpf_set_prec_raw(mpfbuff[mpfstkp], prec);
  else { /* we must increase variable precision */
    mpf_set_prec_raw (mpfbuff[mpfstkp], mpfprec[mpfstkp]);

#if MPTEMP_RENEW
    mpf_clear(mpfbuff[mpfstkp]);
    mpf_init2(mpfbuff[mpfstkp], prec);
#else /* raise precision */
    mpf_set_prec(mpfbuff[mpfstkp], prec);
#endif /* MPTEMP_RENEW */

    mpfprec[mpfstkp] = mpf_get_prec(mpfbuff[mpfstkp]);
  }

  return (tmpf_t) & mpfbuff[mpfstkp++];
}

int
freetmpf(const tmpf_t fp)
{
#if MPTEMP_CHECK
  long prec = mpf_get_prec(fp);
#endif

  if (mpfstkp-- == 0) {
    fprintf(stderr, "Too many deallocations of temporary variables\n");
    exit(EXIT_FAILURE);
  }

#if MPTEMP_DEBUG
  if (fp != mpfbuff[mpfstkp]) {
    fprintf(stderr, "Incorrect deallocation of temporary variable\n");
    exit(EXIT_FAILURE);
  }
#endif

#if MPTEMP_CHECK
  if (prec > mpfprec[mpfstkp]) {
    fprintf(stderr, "Incorrect usage of temporary variable\n");
    exit(EXIT_FAILURE);
  }
#endif /* MPTEMP_CHECK */

  return mpfstkp;
}

void
mpftemp_clear(void)
{
  int i;

#if MPTEMP_DEBUG
  if (mpfstkp != 0) {
    fprintf(stderr, "Some temporary variables are still in use\n");
    exit(EXIT_FAILURE);
  }
#endif

  for (i = 0; i < mpfallp; i++) {
    mpf_set_prec_raw (mpfbuff[i], mpfprec[i]);
    mpf_clear(mpfbuff[i]);
  }

  if (mpfbuff) {
    free(mpfbuff);
    free(mpfprec);
  }

  mpfbuff = NULL;
  mpfprec = NULL;
  mpfsize = 0;
  mpfallp = 0;
  mpfstkp = 0;
}

/***********************************************************
**              mpc_t temporaries                         **
***********************************************************/

static mpc_t *mpcbuff = NULL;
static long  *mpcprec = NULL;
static int mpcsize = 0;
static int mpcstkp = 0;
static int mpcallp = 0;

tmpc_t
gettmpc(long prec)
{
  if (mpcstkp == mpcallp) {
    if (mpcallp == mpcsize) {
      if (mpcsize) {
	mpcbuff = (mpc_t *) realloc(mpcbuff, 2 * mpcsize * sizeof(mpc_t));
	mpcprec = (long *) realloc(mpcprec, 2 * mpcsize * sizeof(long));
	mpcsize *= 2;
      } else {
	mpcbuff = (mpc_t *) malloc(MINBUFFSIZE * sizeof(mpc_t));
	mpcprec = (long *) malloc(MINBUFFSIZE * sizeof(long));
	mpcsize = MINBUFFSIZE;
      }
    }
    mpc_init2(mpcbuff[mpcallp], prec);
    prec = mpc_get_prec(mpcbuff[mpcallp]);
    mpcprec[mpcallp] = prec;
    mpcallp++;
  }

  if (mpcprec[mpcstkp] >= prec)
    mpc_set_prec_raw(mpcbuff[mpcstkp], prec);
  else { /* we must increase variable precision */
    mpc_set_prec_raw (mpcbuff[mpcstkp], mpcprec[mpcstkp]);

#if MPTEMP_RENEW
    mpc_clear(mpcbuff[mpcstkp]);
    mpc_init2(mpcbuff[mpcstkp], prec);
#else /* raise precision */
    mpc_set_prec(mpcbuff[mpcstkp], prec);
#endif /* MPTEMP_RENEW */

    mpcprec[mpcstkp] = mpc_get_prec(mpcbuff[mpcstkp]);
  } /* precision check */

  return (tmpc_t) & mpcbuff[mpcstkp++];
}

int
freetmpc(const tmpc_t fp)
{
#if MPTEMP_CHECK
  long prec = mpc_get_prec(fp);
#endif

  if (mpcstkp-- == 0) {
    fprintf(stderr, "Too many deallocations of temporary variables\n");
    exit(EXIT_FAILURE);
  }

#if MPTEMP_DEBUG
  if (fp != mpcbuff[mpcstkp]) {
    fprintf(stderr, "Incorrect deallocation of temporary variable\n");
    exit(EXIT_FAILURE);
  }
#endif

#if MPTEMP_CHECK
  if (prec > mpcprec[mpcstkp]) {
    fprintf(stderr, "Incorrect usage of temporary variable\n");
    exit(EXIT_FAILURE);
  }
#endif /* MPTEMP_CHECK */

  return mpcstkp;
}

void
mpctemp_clear(void)
{
  int i;

#if MPTEMP_DEBUG
  if (mpcstkp != 0) {
    fprintf(stderr, "Some temporary variable are still in use\n");
    exit(EXIT_FAILURE);
  }
#endif

  for (i = 0; i < mpcallp; i++) {
    mpc_set_prec_raw (mpcbuff[i], mpcprec[i]);
    mpc_clear(mpcbuff[i]);
  }

  if (mpcbuff) {
    free(mpcbuff);
    free(mpcprec);
  }

  mpcbuff = NULL;
  mpcprec = NULL;
  mpcsize = 0;
  mpcallp = 0;
  mpcstkp = 0;
}

/***********************************************************
**              global functions                          **
***********************************************************/

void
mptemp_clear(void)
{
  mpftemp_clear();
  mpctemp_clear();
}

#endif /* NOMPTEMP */
