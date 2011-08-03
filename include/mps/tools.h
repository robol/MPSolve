/***********************************************************
**                        C Tools                         **
**                      Version 2.0                       **
**                                                        **
**             Written by Giuseppe Fiorentino             **
**                 (fiorent@dm.unipi.it)                  **
***********************************************************/

#ifndef MPS_TOOLS_H_
#define MPS_TOOLS_H_

#ifdef __cplusplus
extern "C"
{
#endif

/* consts */
#define LOG2     0.69314718055994530941
#define LOG10    2.30258509299404568401
#define LOG10_2  0.30102999566398119521
#define LOG2_10  3.32192809488736234787
#define PI       3.14159265358979323846

/* macros */
#define MAX(A, B)  ( (A) > (B) ? (A) : (B) )
#define MIN(A, B)  ( (A) < (B) ? (A) : (B) )

/* types */
#ifndef __USE_BOOL_AS_BOOLEAN
  typedef enum
  { false = 0, true = 1 } mps_boolean;
#else
  /* Small workaround to make matlab module work; there is,
   * int matlab headers, already a false keyword defined, so
   * reusing it here make compilation fail. */
  typedef bool mps_boolean;
#endif				/* mps_boolean */

#define mps_boolean_to_string(x) ((x) == true) ? "true" : "false"

/* functions */
  void randomize (unsigned int seed);
  double drand (void);
  double dbl_set_2dl (double d, long int l);
  void dbl_get_2dl (double *rd, long int *rl, double d);
  double dbl_get_mant (double d);
  int dbl_get_exp (double d);

/* vector support functions */
#define mps_boolean_valloc(N)		(mps_boolean *) malloc((N) * sizeof(mps_boolean))
  void mps_boolean_vinit (mps_boolean v[], unsigned long int size);
#define mps_boolean_vclear(V, N)		mps_boolean_vinit(V, N)
#define mps_boolean_vfree(V)		free(V)

/* vector support functions */
#define char_valloc(N)			(char *) malloc((N) * sizeof(char))
  void char_vinit (char v[], unsigned long int size);
#define char_vclear(V, N)		char_vinit(V, N)
#define char_vfree(V)			free(V)

#define int_valloc(N)			(int *) malloc((N) * sizeof(int))
  void int_vinit (int v[], unsigned long int size);
#define int_vclear(V, N)		int_vinit(V, N)
#define int_vfree(V)			free(V)

#define long_valloc(N)			(long *) malloc((N) * sizeof(long))
  void long_vinit (long v[], unsigned long int size);
#define long_vclear(V, N)		lng_vinit(V, N)
#define long_vfree(V)			free(V)

#define float_valloc(N)			(float *) malloc((N) * sizeof(float))
  void float_vinit (float v[], unsigned long int size);
#define float_vclear(V, N)		float_vinit(V, N)
#define float_vfree(V)			free(V)

#define double_valloc(N)		(double *) malloc((N) * sizeof(double))
  void double_vinit (double v[], unsigned long int size);
#define double_vclear(V, N)		double_vinit(V, N)
#define double_vfree(V)			free(V)

/*
 * End of extern "C" {
 *   ...
 * }
 */
#ifdef __cplusplus
}
#endif

#endif
