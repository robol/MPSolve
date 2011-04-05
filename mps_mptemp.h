/***********************************************************
**        Memory Manager for mpf and mpc Temporaries      **
**                      Version 1.0                       **
**                                                        **
**             Written by Giuseppe Fiorentino             **
**                 (fiorent@dm.unipi.it)                  **
***********************************************************/

#ifndef __MPTEMP_H__
#define __MPTEMP_H__

/* include GMP header file */
#include "gmp.h"
#include "mps_mpc.h"

#ifndef NOMPTEMP  /* enable MPTEMP variables              */

/***********************************************************
**                    MPTEMP ENABLED                      **
***********************************************************/

/* tmpf_t function prototypes */

typedef __mpf_struct *tmpf_t;

#define mpftemp_init()            	mpftemp_clear()
void mpftemp_clear(void);

tmpf_t gettmpf(long int prec);
int freetmpf(const tmpf_t fp);

#define tmpf_init2(F, P)          	F = gettmpf(P)
#define tmpf_init3(F, P, I)       	F = gettmpf2(P, &I)
void tmpf_raise_prec(tmpf_t f, long int prec, int pos);
#define tmpf_clear(F)             	freetmpf(F)

/* tmpc_t function prototypes */

typedef __mpc_struct *tmpc_t;

#define mpctemp_init()            	mpctemp_clear()
void mpctemp_clear(void);

tmpc_t gettmpc(long int prec);
int freetmpc(const tmpc_t cp);

#define tmpc_init2(C, P)          	C = gettmpc(P)
#define tmpc_clear(C)             	freetmpc(C)

/* global function prototypes */

#define mptemp_init()             	mptemp_clear()
void mptemp_clear(void);

/***********************************************************
**                                                        **
***********************************************************/

#else /* NOMPTEMP defined - disable MPTEMP variables      */

/***********************************************************
**                   MPTEMP DISABLED                      **
***********************************************************/

/* simulated tmpf_t function prototypes */

typedef mpf_t tmpf_t;

#define mpftemp_init()			/* void */
#define mpftemp_clear()			/* void */

#define tmpf_init2(F, P)		mpf_init2(F, P)
#define tmpf_init3(F, P, I)		mpf_init2(F, P)
#define tmpf_raise_prec(F, P, I) 	mpf_set_prec(F, P)
#define tmpf_clear(F)			mpf_clear(F)

/* tmpc_t function prototypes */

typedef mpc_t tmpc_t;

#define mpctemp_init()			/* void */
#define mpctemp_clear()			/* void */

#define tmpc_init2(C, P)		mpc_init2(C, P)
#define tmpc_clear(C)			mpc_clear(C)

/* global function prototypes */

#define mptemp_init()			/* void */
#define mptemp_clear()			/* void */

#endif /* ndef NOMPTEMP */

#endif /* def __MPTEMP_H__ */
