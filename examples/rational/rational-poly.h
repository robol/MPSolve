#include <mps/mps.h>

#define MPS_RATIONAL_POLY(t) (MPS_POLYNOMIAL_CAST (mps_rational_poly, t))
#define MPS_IS_RATIONAL_POLY(t) (mps_polynomial_check_type (t, "mps_rational_poly"))

struct mps_rational_poly {
  mps_polynomial p;

  /* Exponent of the denominators. */
  int pp;

  /* Number of summands in the rational function. */
  int nn;

  /* Pointer to the coefficients of the summands. */
  cplx_t * a; 
  cplx_t * b;
  cplx_t * c;

  /* Leading coefficient of the polynomial */
  mpcf_t lc;
};

typedef struct mps_rational_poly mps_rational_poly;

mps_rational_poly * mps_rational_poly_new (mps_context * ctx, cplx_t * a,
					   cplx_t * b, cplx_t * c, int nn,
					   int pp);

mps_boolean mps_rational_poly_feval (mps_context * ctx, mps_polynomial * p, 
				     cplx_t x, cplx_t value, double * error);
mps_boolean mps_rational_poly_deval (mps_context * ctx, mps_polynomial * p, 
				     cdpe_t x, cdpe_t value, rdpe_t error);
mps_boolean mps_rational_poly_meval (mps_context * ctx, mps_polynomial * p, 
				     mpcf_t x, mpcf_t value, rdpe_t error);

void mps_rational_poly_fnewton (mps_context * ctx, mps_polynomial * p,
				mps_approximation * root, cplx_t corr);
void mps_ratioanl_poly_dnewton (mps_context * ctx, mps_polynomial * p,
				mps_approximation * root, cdpe_t corr);
void mps_rational_poly_mnewton (mps_context * ctx, mps_polynomial * p,
				mps_approximation * root, mpcf_t corr, long int wp);

void mps_rational_poly_free (mps_context * ctx, mps_rational_poly * p);

void mps_rational_poly_get_leading_coefficient (mps_context * ctx, 
						mps_polynomial * p, 
						mpcf_t lc);
