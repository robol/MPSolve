#include "quadratic-poly.h"

/** 
 * @brief Allocate a new Quadratic polynomial of the specified level. 
 * 
 * Quadratic polynomials are defined by recurrence as \f$P_0(z) = 0\f$, 
 * \f$P_{k+1} = P_{k}(z)^2 + z\f$. The polynomial \f$P_k(z)\f$ is called the
 * \f$k\f$-th level Quadratic polynomial. 
 *
 * @param ctx The current mps_context. 
 * @param level The desired level for the Quadratic polynomial. 
 * @return A pointer to a newly-allocated mps_quadratic_poly. 
 */
mps_quadratic_poly *
mps_quadratic_poly_new (mps_context * ctx, int level)
{
  mps_quadratic_poly *mp = mps_new (mps_quadratic_poly);
  mps_polynomial *p = MPS_POLYNOMIAL (mp);
  mps_polynomial_init (ctx, p);

  mp->level = level;
  p->degree = pow(2, level) - 1;

  p->type_name = "mps_quadratic_poly";
  p->structure = MPS_STRUCTURE_REAL_INTEGER;
  p->density = MPS_DENSITY_DENSE;
  
  /* Load methods */
  p->fnewton = mps_quadratic_poly_fnewton;
  p->dnewton = mps_quadratic_poly_dnewton;
  p->mnewton = mps_quadratic_poly_mnewton;
  p->feval = mps_quadratic_poly_feval;
  p->deval = mps_quadratic_poly_deval;
  p->meval = mps_quadratic_poly_meval;
  p->dstart = mps_quadratic_poly_dstart;
  p->fstart = mps_quadratic_poly_fstart;
  p->free = mps_quadratic_poly_free;

  return mp;
}

/**
 * @brief Free a Quadratic polynomial. 
 *
 * @param ctx The current mps_context. 
 * @param p The polynomial that has to be freed. 
 */
void
mps_quadratic_poly_free (mps_context * ctx, 
			  mps_polynomial * p)
{
  /* Quadratic polynomials haven't any data, so there's essentially
   * nothing to do here. */
  free (p);
}
