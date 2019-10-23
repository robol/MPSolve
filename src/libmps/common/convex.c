/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@sns.it>
 */


/*
 * Routines for the computation of convex hulls.
 *
 * More specifically, the routine <code>fconvex</code> is the one
 * that is intentended to be used and, given a vector of double
 * <code>a</code>, computes the convex hull of the set
 * \f[ \mathcal{S} = \{ (i, a_i) \ | \ i \in \{ 1, \dots, n\} \} \f]
 *
 * The computed set is returned in a vector of mps_booleans <code>h</code>
 * such that its vertices are
 * \f[ \mathcal{V} = \{ (i, a_i) \ | \ i \in \{ 1, \dots, n \} \ \text{and} \ h_i \ \text{is true} \ \} \f]
 */

#include <mps/mps.h>
#include <string.h>

static const double TOLER = 0.4;       /* slope tolerace */

/**
 * @brief find max lo<j<=i : h[j]
 *
 * More clearly, find the minimum index \f$j\f$ such that
 * \f$lo < j \leq i\f$ and \f$(j, a_j)\f$ is a vertex of
 * the convex hull of the points
 * \f[ \{ (k, a_k) \ | \ k \in \{ lo, \dots, i \} \} \f]
 */
static int
mps_left (mps_context * s, int i, int lo, int * h)
{
  if (i == lo)
    return lo;
  for (i--; i > lo; i--)
    if (h[i])
      break;
  return i;
}

/**
 * @brief find min i<=j<up : h[j]
 *
 * More clearly, find the maximum index \f$j\f$ such that
 * \f$i \leq j < up\f$ and \f$(j, a_j)\f$ is a vertex of
 * the convex hull of the points
 * \f[ \{ (k, a_k) \ | \ k \in \{ i, \dots, up \} \} \f]
 */
int
mps_right (mps_context * s, int i, int up, int * h)
{
  if (i == up)
    return up;
  for (i++; i < up; i++)
    if (h[i])
      break;
  return i;
}

/**
 * @brief convexity test of the points \f$ \{ (il, a_{il}), (i, a_i), (ir, a_{ir})\}\f$.
 *
 * Check if the points in the given set are "enough convex", i.e. if
 * the set is \f$ \{ P_1, P_2, P_3 \} \f$ check if the slope
 * of the line from \f$P_1\f$ to \f$P_2\f$ is at least <code>TOLER</code>
 * less than the slope of the line joining \f$P_2\f$ and \f$P_3\f$.
 *
 * @param s The <code>mps_context</code> associated with the current computation.
 * @param il index of the first point
 * @param i  index of the middle point
 * @param ir index of the last point
 * @param a  array with the points
 */
static mps_boolean
mps_fctest (mps_context * s, int il, int i, int ir, double a[])
{
  double s1, s2;

  s1 = (a[i] - a[il]) * (ir - i);
  s2 = (a[ir] - a[i]) * (i - il);
  /*#rimosso giugno 2000
     return (s1 - s2 > (i-il)*(ir-i)*TOLER);
   */
  return(s1 - s2 > TOLER);
}

/**
 * @brief Merge two adjacent convex hulls [lo, i] and [i, hi].
 *
 * @param s The <code>mps_context</code> associated with the current computation.
 * @param lo starting index of the points of the first convex hull
 * in the vector <code>a</code>.
 * @param i last index of the points in the first convex hull and
 * first of index of the points in the second convex hull.
 * @param up last index of the points in the second convex hull
 * @param a array of points
 */
static void
mps_fmerge (mps_context * s, int lo, int i, int up, double a[], int * h)
{
  int il, ir, ill, irr;
  mps_boolean tstl, tstr;

  ill = lo;
  irr = up;
  il = mps_left (s, i, lo, h);
  ir = mps_right (s, i, up, h);
  if (mps_fctest (s, il, i, ir, a))
    return;
  h[i] = false;
  do
    {
      if (il == lo)
        tstl = true;
      else
        {
          ill = mps_left (s, il, lo, h);
          tstl = mps_fctest (s, ill, il, ir, a);
        }
      if (ir == up)
        tstr = true;
      else
        {
          irr = mps_right (s, ir, up, h);
          tstr = mps_fctest (s, il, ir, irr, a);
        }
      if (!tstl)
        {
          h[il] = false;
          il = ill;
        }
      if (!tstr)
        {
          h[ir] = false;
          ir = irr;
        }
    } while (!(tstl && tstr));
}

/**
 * @brief compute the convex hull of the data set a[].
 *
 * The result
 * is in the mps_boolean vector <code>h[]</code>. The algorithm successively
 * merges adjacent convex hulls of sizes 2, 4, 8, ...
 *
 * @param s The <code>mps_context</code> associated with the current computation.
 * @param a vector of points whose convex hull must be computed.
 * @param n size of the vector <code>a</code>.
 *
 * @return h An integer array that is 1 only in the positions corresponding to
 * selected vertexes. 
 */
MPS_PRIVATE int *
mps_fconvex (mps_context * s, int n, double a[])
{
  int m, c;
  
  int * h = mps_newv (int, (s->n + 1));
  memset (h, 1, sizeof (int) * (s->n + 1));

  for (m = 1; m < s->n; m <<= 1)
    for (c = m; c < s->n; c += 2 * m)
      mps_fmerge (s, c - m, c, MIN (s->n, c + m), a, h);

  return h;
}

MPS_PRIVATE mps_vertex *
mps_vertex_copy (mps_context * ctx, mps_vertex * source)
{
  mps_vertex * copy = mps_new (mps_vertex);
  memmove (copy, source, sizeof (mps_vertex));
  return copy;
}


MPS_PRIVATE mps_linear_hypograph * 
mps_linear_hypograph_new (mps_context * ctx)
{
  mps_linear_hypograph * sl =  mps_new (mps_linear_hypograph);

  sl->first = sl->last = NULL;
  sl->n = 0;

  return sl;
}

MPS_PRIVATE void
mps_linear_hypograph_free (mps_context * ctx, mps_linear_hypograph * sl)
{
  free (sl);
}

MPS_PRIVATE void 
mps_linear_hypograph_append_node (mps_context * ctx, mps_linear_hypograph * sl, 
				 mps_vertex * node)
{
  if (sl->last)
    {
      node->previous = sl->last;
      sl->last->next = node;
      sl->last = node;
    }
  else
    {
      node->previous = NULL;
      sl->last = sl->first = node;
    }

  node->next = NULL;
}

MPS_PRIVATE void
mps_linear_hypograph_remove_node (mps_context * ctx, mps_linear_hypograph * sl, 
				 mps_vertex * node)
{
  mps_vertex * previous = node->previous;
  mps_vertex * next = node->next;

  if (previous)
    previous->next = next;

  if (next)
    next->previous = previous;

  if (node == sl->first)
    sl->first = next;

  if (node == sl->last)
    sl->last = previous;
}

MPS_PRIVATE mps_linear_hypograph * 
mps_linear_hypograph_concat (mps_context * ctx, mps_linear_hypograph * sl1,
			    mps_linear_hypograph * sl2)
{
  sl1->last->next = sl2->first;
  sl1->last = sl2->last;
  sl2->first->previous = sl1->last;

  return sl1;
}

MPS_PRIVATE mps_linear_hypograph *
mps_linear_hypograph_detach_tail (mps_context * ctx, mps_linear_hypograph * sl, 
				 mps_vertex * node)
{
  mps_linear_hypograph * sl_tail = mps_linear_hypograph_new (ctx);
  mps_vertex * new_first = mps_vertex_copy (ctx, node);

  sl_tail->last = sl->last;
  sl_tail->first = new_first;
  new_first->previous = NULL;
  
  node->next = NULL;
  sl->last = node;

  return sl_tail;
}

MPS_PRIVATE mps_linear_hypograph *
mps_convex_hull (mps_context * s, mps_linear_hypograph * l)
{
  return NULL;
}


