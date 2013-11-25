#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

/* Verify that a cluster can be correctly created 
 * and the number of roots in it is correct even after
 * adding other roots in it. */
START_TEST (cluster_create)
{
  mps_context *s = mps_context_new ();
  mps_cluster *cluster = mps_cluster_empty (s);

  // Add some roots to the cluster and verify that
  // the correct number of roots is maintained over time. 
  fail_unless (cluster->n == 0, "An empty cluster should have 0 roots");

  mps_cluster_insert_root (s, cluster, 45);
  fail_unless (cluster->n == 1, "A cluster with a root should have"
               " cluster->n == 1");

  mps_cluster_free (s, cluster);
  mps_context_free (s);
}
END_TEST

/* Verify proper handling of isolation in cluster. 
 * A fake cluster set is created for the polynomial
 * (x - 2)^2 * (x - 1) with the approximations: 
 * - 1.01
 * - 2.01
 * - 1.99
 * The clusters obtained from these approximations should
 * be two: one containing the approximation 1.01 and the
 * other two containing 2.01 and 1.99. 
 */
START_TEST (cluster_isolation)
{
  mps_cluster_item *cluster_item = NULL;
  mps_cluster *cluster = NULL;
  mps_context *s = mps_context_new ();
  
  // Set the input polynomial that we have chosen, i.e.
  // x^3 - 5x^2 + 8x - 4
  mps_monomial_poly *p = mps_monomial_poly_new (s, 3);
  mps_monomial_poly_set_coefficient_int (s, p, 3, 1, 0);
  mps_monomial_poly_set_coefficient_int (s, p, 2, -5, 0);
  mps_monomial_poly_set_coefficient_int (s, p, 1, 8, 0);
  mps_monomial_poly_set_coefficient_int (s, p, 0, -4, 0);

  mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));
  mps_allocate_data (s);

  // Select the starting approximations
  cplx_set_d (s->root[0]->fvalue, 2.01, 0.0);
  cplx_set_d (s->root[1]->fvalue, 1.99, 0.0);
  cplx_set_d (s->root[2]->fvalue, 1.01, 0.0);

  // Set the Newton radii to a sufficient big value such that 
  // it won't disturb the Gerschgorin ones. 
  s->root[0]->frad = s->root[1]->frad = s->root[2]->frad = 5.0;

  // Compute the gerschgorin radii of inclusions
  double * gerschgorin_radii = mps_newv (double, 3);
  mps_fradii (s, s->active_poly, gerschgorin_radii);

  mps_fcluster (s, gerschgorin_radii, 2.0 * s->n);

  // Check that we have two clusters and that are the
  // cluster that we are expecting. 
  fail_unless (s->clusterization->n == 2, "There should be two clusters in"
               " the given example, but %d were found", s->clusterization->n);

  for (cluster_item = s->clusterization->first; cluster_item; 
       cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;
      if (cluster->n == 1)
        fail_unless (cplx_Re (s->root[cluster->first->k]->fvalue) == 1.01, 
                     "The isolated approximation in the example should be 1.01");
      else
        {
          double first_real_part = cplx_Re (s->root[cluster->first->k]->fvalue);
          double second_real_part = 
            cplx_Re (s->root[cluster->first->next->k]->fvalue);
          fail_unless ((first_real_part == 2.01 && second_real_part == 1.99) ||
                       (first_real_part == 1.99 && second_real_part == 2.01),
                       "The approximations in the cluster with cardinality two"
                       " should be 1.99 and 2.01");
        }
    }

  free (gerschgorin_radii);
  mps_context_free (s);
}
END_TEST


int 
main (void) 
{
  int number_failed;
  starting_setup ();

  Suite *s = suite_create ("Cluster analsysis");
  TCase *tc_management = tcase_create ("Cluster management");

  // Add tests of the Cluster management test case
  tcase_add_test (tc_management, cluster_create);
  tcase_add_test (tc_management, cluster_isolation);

  suite_add_tcase (s, tc_management);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
