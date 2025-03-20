#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

START_TEST (test_version_strings)
{
  char buffer[255];

  /* We check that all the different ways of producing the
   * string describing MPSolve's version actually match */
  sprintf(buffer, "%u.%u.%u", mps_get_major_version(),
      mps_get_minor_version(), mps_get_patch_version());

  ck_assert (strcmp(buffer, mps_get_version()) == 0);

  sprintf(buffer, "%u.%u.%u", MPS_MAJOR_VERSION,
      MPS_MINOR_VERSION, MPS_PATCH_VERSION);

  ck_assert (strcmp(buffer, mps_get_version()) == 0);
}
END_TEST

int main()
{
  int number_failed = 0;

  starting_setup ();

  Suite *s = suite_create ("Version");
  TCase *tc_version = tcase_create ("Version number");

  tcase_add_test (tc_version, test_version_strings);

  suite_add_tcase (s, tc_version);

  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return(number_failed != 0);
}
