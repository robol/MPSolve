/***********************************************************
**       Multiprecision Polynomial Solver (MPSolve)       **
**                 Version 2.2, May 2001                  **
**                                                        **
**                      Written by                        **
**       Dario Andrea Bini and Giuseppe Fiorentino        **
**       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
**                                                        **
** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

/**
 * @file
 * @brief Routines for the computation of convex hulls.
 */

#include "mps.h"

const double TOLER = 0.4;	/* slope tolerace */

/***********************************************************
 *                      SUBROUTINE LEFT                    *
 ***********************************************************
   find max lo<j<=i : h[j]
 **********************************************************/
int
left(int i, int lo)
{
  if (i == lo)
    return lo;
  for (i--; i > lo; i--)
    if (h[i])
      break;
  return i;
}

/***********************************************************
 *                      SUBROUTINE RIGHT                   *
 ***********************************************************
   find min i<=j<=up : h[j]
 **********************************************************/
int
right(int i, int up)
{
  if (i == up)
    return up;
  for (i++; i < up; i++)
    if (h[i])
      break;
  return i;
}

/***********************************************************
 *                      SUBROUTINE FCTEST                  *
 ***********************************************************
   convexity test {(il, a[il]), (i, a[i]), (ir, a[ir])}
 **********************************************************/
boolean
fctest(int il, int i, int ir, double a[])
{
  double s1, s2;

  s1 = (a[i] - a[il]) * (ir - i);
  s2 = (a[ir] - a[i]) * (i - il);
  /*#rimosso giugno 2000
  return (s1 - s2 > (i-il)*(ir-i)*TOLER);
  */
  return (s1 - s2 > TOLER);
}

/**
 * @brief Merge two adjacent convex hulls [lo, i] and [i, hi].
 * @param lo starting index of the points of the first convex hull
 * in the vector <code>a</code>.
 * @param i last index of the points in the first convex hull and
 * first of index of the points in the second convex hull.
 * @param up last index of the points in the second convex hull
 * @param a array of points
 */
void
fmerge(int lo, int i, int up, double a[])
{
  int il, ir, ill, irr;
  boolean tstl, tstr;

  ill = lo;
  irr = up;
  il = left(i, lo);
  ir = right(i, up);
  if (fctest(il, i, ir, a))
    return;
  h[i] = false;
  do {
    if (il == lo)
      tstl = true;
    else {
      ill = left(il, lo);
      tstl = fctest(ill, il, ir, a);
    }
    if (ir == up)
      tstr = true;
    else {
      irr = right(ir, up);
      tstr = fctest(il, ir, irr, a);
    }
    if (!tstl) {
      h[il] = false;
      il = ill;
    }
    if (!tstr) {
      h[ir] = false;
      ir = irr;
    }
  }
  while (!(tstl && tstr));
}

/**
 * @brief compute the convex hull of the data set a[].
 *
 * The result
 * is in the boolean vector <code>h[]</code>. The algorithm successively
 * merges adjacent convex hulls of sizes 2, 4, 8, ...
 *
 * @param a vector of points whose convex hull must be computed.
 * @param n size of the vector <code>a</code>.
 */
void
fconvex(int n, double a[])
{
  int m, c;

  for (m = 0; m <= n; m++)
    h[m] = true;

  for (m = 1; m < n; m <<= 1)
    for (c = m; c < n; c += 2 * m)
      fmerge(c - m, c, MIN(n, c + m), a);
}
