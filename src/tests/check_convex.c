#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

/* Check that we can correctly create and destroy a hypograph. */
START_TEST (hypograph_creation)
{
  mps_context * ctx = mps_context_new ();
  mps_linear_hypograph * l = mps_linear_hypograph_new (ctx);

  mps_linear_hypograph_free (ctx, l);
  mps_context_free (ctx);
}
END_TEST

int
main (void)
{
  int number_failed;

  starting_setup ();

  Suite *s = suite_create ("Convex operations");
  TCase *tc_hypograph = tcase_create ("Hypograph management");

  // Hypograph management
  suite_add_tcase (s, tc_hypograph);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return(number_failed != 0);
}
