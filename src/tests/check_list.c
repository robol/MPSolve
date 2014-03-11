#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

START_TEST (list_creation)
{
  mps_list * list = mps_list_new ();

  fail_unless (mps_list_size (list) == 0, 
	       "Newly created list's size differs from 0.");

  mps_list_free (list);
}
END_TEST

START_TEST (list_foreach)
{
  mps_list * list = mps_list_new ();
  int test_elements[] = { 5,2,12,2,12 };
  int i, n = 5;
  int * ptr;

  for (i = 0; i < n; i++)
    mps_list_append (list, mps_list_element_new (test_elements + i));
  
  i = 0;
  MPS_LIST_FOREACH (int, ptr, list)
    {
      fail_unless (*ptr == test_elements[i++], 
		   "List elements differs from the array that was used to fill it.");
    }
  
}
END_TEST

START_TEST (element_creation)
{
  mps_list_element * el = mps_list_element_new (NULL);
  mps_list_element_free (el);
}
END_TEST

int
main (void)
{
  int number_failed;

  starting_setup ();

  Suite *s = suite_create ("Containers");
  TCase *tc_elements = tcase_create ("Basic list elements operations");
  TCase *tc_lists = tcase_create ("List operations");

  // Basic operations on list elements
  tcase_add_test (tc_elements, element_creation);

  // Basic operations on lists
  tcase_add_test (tc_lists, list_creation);
  tcase_add_test (tc_lists, list_foreach);

  suite_add_tcase (s, tc_elements);
  suite_add_tcase (s, tc_lists);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return(number_failed != 0);
}
