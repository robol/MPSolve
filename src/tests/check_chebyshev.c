#include <check.h>
#include <mps/mps.h>

START_TEST (test_chebyshev_poly_80)
{
  mps_context * ctx = mps_context_new ();
  mps_chebyshev_poly *cp = mps_chebyshev_poly_new (ctx, 80, MPS_STRUCTURE_REAL_INTEGER);
  mpc_t *mroots = NULL;
  rdpe_t *radii = NULL;
  mpq_t one, zero;
  int i;

  mpq_init (one);
  mpq_init (zero);

  mpq_set_ui (one, 1U, 1U);
  mpq_set_ui (zero, 0U, 1U);

  mps_chebyshev_poly_set_coefficient_q (ctx, cp, 80, one, zero);
  for (i = 0; i < 80; i++)
    {
      mps_chebyshev_poly_set_coefficient_q (ctx, cp, i, zero, zero);
    }

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (cp));
  mps_context_select_algorithm (ctx, MPS_ALGORITHM_SECULAR_GA);
  mps_thread_pool_set_concurrency_limit (ctx, NULL, 1);
  mps_mpsolve (ctx);

  /* Check that the roots are correct */
  mps_context_get_roots_m (ctx, &mroots, &radii);
  for (i = 0; i < 80; i++)
    {
      int j, found_root = -1; 
      rdpe_t expected_root;
      rdpe_t diff;
      double epsilon = DBL_MAX;

      rdpe_set_d (expected_root, cos ( (2.0*i+1 ) / 160 * PI ));; 

      printf ("[Chebyshev tests] Expected root %d: (%1.20lf, 0)\n",
	      i, rdpe_get_d (expected_root));

      for (j = 0; j < 80; j++)
	{
	  cdpe_t ctmp; 
	  double residue;

	  mpc_get_cdpe (ctmp, mroots[j]);
	  rdpe_sub (diff, expected_root, cdpe_Re (ctmp));

	  residue = sqrt( pow (fabs (rdpe_get_d (diff)), 2) + pow (fabs( rdpe_get_d (cdpe_Im (ctmp))), 2) );
	  if (residue < epsilon)
	    {
	      epsilon = residue;
	      found_root = j;
	    }
	}

      printf ("[Chebyshev tests] Residue for approximation %3d: %e\n", i, epsilon);
      printf ("[Chebyshev tests] Inclusion radii: %e\n", rdpe_get_d (radii[found_root]));
      fail_unless (epsilon < 4.0 * DBL_EPSILON + rdpe_get_d (radii[found_root]));
    }

  mps_polynomial_free (ctx, MPS_POLYNOMIAL (cp));
  mps_context_free (ctx);

  mpq_clear (one);
  mpq_clear (zero);
}
END_TEST

START_TEST (test_chebyshev_poly_20)
{
  mps_context * ctx = mps_context_new ();
  mps_chebyshev_poly *cp = mps_chebyshev_poly_new (ctx, 20, MPS_STRUCTURE_REAL_INTEGER);
  mpc_t *mroots = NULL;
  rdpe_t *radii = NULL;
  mpq_t one, zero;
  int i;

  mpq_init (one);
  mpq_init (zero);

  mpq_set_ui (one, 1U, 1U);
  mpq_set_ui (zero, 0U, 1U);

  mps_chebyshev_poly_set_coefficient_q (ctx, cp, 20, one, zero);
  for (i = 0; i < 20; i++)
    {
      mps_chebyshev_poly_set_coefficient_q (ctx, cp, i, zero, zero);
    }

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (cp));
  mps_context_select_algorithm (ctx, MPS_ALGORITHM_SECULAR_GA);
  mps_thread_pool_set_concurrency_limit (ctx, NULL, 1);
  mps_mpsolve (ctx);

  /* Check that the roots are correct */
  mps_context_get_roots_m (ctx, &mroots, &radii);
  for (i = 0; i < 20; i++)
    {
      int j, found_root = -1; 
      rdpe_t expected_root;
      rdpe_t diff;
      double epsilon = DBL_MAX;

      rdpe_set_d (expected_root, cos ( (2.0*i+1 ) / 40 * PI ));; 

      printf ("[Chebyshev tests] Expected root %d: (%1.20lf, 0)\n",
	      i, rdpe_get_d (expected_root));

      for (j = 0; j < 20; j++)
	{
	  cdpe_t ctmp; 
	  double residue;

	  mpc_get_cdpe (ctmp, mroots[j]);
	  rdpe_sub (diff, expected_root, cdpe_Re (ctmp));

	  residue = sqrt( pow (fabs (rdpe_get_d (diff)), 2) + pow (fabs( rdpe_get_d (cdpe_Im (ctmp))), 2) );
	  if (residue < epsilon)
	    {
	      epsilon = residue;
	      found_root = j;
	    }
	}

      printf ("[Chebyshev tests] Residue for approximation %3d: %e\n", i, epsilon);
      printf ("[Chebyshev tests] Inclusion radii: %e\n", rdpe_get_d (radii[found_root]));
      fail_unless (epsilon < 4.0 * DBL_EPSILON + rdpe_get_d (radii[found_root]));
    }

  mps_polynomial_free (ctx, MPS_POLYNOMIAL (cp));
  mps_context_free (ctx);

  mpq_clear (one);
  mpq_clear (zero);
}
END_TEST

Suite*
chebyshev_suite (void)
{
  Suite *s = suite_create ("chebyshev");

  TCase *tcase_t = tcase_create ("Solution of Chebyshev polynomials");
  tcase_add_test (tcase_t, test_chebyshev_poly_20);
  tcase_add_test (tcase_t, test_chebyshev_poly_80);

  suite_add_tcase (s, tcase_t);
  return s;
}

int 
main (void)
{
  Suite *cs = chebyshev_suite ();
  SRunner *sr = srunner_create (cs);
  int number_failed;

  srunner_run_all (sr, CK_NORMAL);

  /* Get number of failed test and report */
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
