#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

START_TEST (basics_addition)
{
  mpcf_t a, b, c;
  cplx_t x;

  mpcf_init2 (a, 150 * LOG2_10 + DBL_MANT_DIG);
  mpcf_init2 (b, 150 * LOG2_10 + DBL_MANT_DIG);
  mpcf_init2 (c, 150 * LOG2_10 + DBL_MANT_DIG);

  mpcf_set_d (a, 1.0, 0.0);
  mpcf_set_d (b, 1e-150, 0.0);

  mpcf_add (c, a, b);
  mpcf_sub_eq (c, a);

  mpcf_get_cplx (x, c);

  fail_unless (fabs (cplx_Re (x) - 1e-150) < DBL_EPSILON * 1e-150 * 8.0,
               "Loss of precision in multiprecision arithmetic: addition");

  mpcf_clear (a);
  mpcf_clear (b);
  mpcf_clear (c);
}
END_TEST

START_TEST (basics_multiplication)
{
  mpcf_t a, b, c;
  long int precisions[] = { 64, 128, 1024 };
  int i;

  for (i = 0; i < 3; i++)
    {
      mpcf_init2 (a, precisions[i]);
      mpcf_init2 (b, precisions[i]);
      mpcf_init2 (c, precisions[i]);

      /* Test some basics operations that may cause overflow in standard floating
       * point. We check if (1 + 2i) * (1 + epsilon) / (1 + epsilon) == (1 + 2i) */
      mpcf_set_ui (a, 1U, 2U);
      mpcf_set_ui (b, 1U, 0U);

      /* Construct epsilon */
      {
        cdpe_t epsilon;

        cdpe_set (epsilon, cdpe_zero);
        rdpe_set_2dl (cdpe_Re (epsilon), 1.0, -precisions[i] + DBL_MANT_DIG);

        mpcf_set_cdpe (c, epsilon);
        mpcf_add_eq (b, c);
      }

      /* Perform the multiplication */
      mpcf_mul (c, a, b);
      mpcf_div (c, c, b);

      /* Check if the result is correct in floating point. That should be true since epsilon
       * was greater than 2.0^{ - precision + DBL_MANT_DIG }. */
      {
        cplx_t x, y;

        mpcf_get_cplx (x, c);
        cplx_set_d (y, 1.0, 2.0);

        cplx_sub_eq (x, y);

        fail_unless (cplx_mod (x) < cplx_mod (y) * 8.0 * DBL_EPSILON,
                     "Floating point errors in complex multiplications");
      }

      mpcf_clear (a);
      mpcf_clear (b);
      mpcf_clear (c);
    }
}
END_TEST

int
main (void)
{
  int number_failed;

  starting_setup ();

  Suite *s = suite_create ("Multiprecision arithmetic");
  TCase *tc_basics = tcase_create ("Simple arithmetic");

  // Basic operations
  tcase_add_test (tc_basics, basics_addition);
  tcase_add_test (tc_basics, basics_multiplication);

  suite_add_tcase (s, tc_basics);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return(number_failed != 0);
}
