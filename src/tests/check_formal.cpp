#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

#include <iostream>

START_TEST (monomial_creation)
{
  mps::formal::Monomial one("1", 0);
  fail_unless (one.coefficientReal() == 1,
	       "Monomial(\"1\", 0) != 1");

  mps::formal::Monomial exp1("1.01e2", 12);
  fail_unless ((exp1.coefficientReal() == 101) && 
	       (exp1.degree() == 12),
	       "Monomial(\"1.01e2\",12) != 101x^12");
	       
}
END_TEST

START_TEST (monomial_floating)
{
  mps::formal::Monomial test("0.25", 0);
  fail_unless (test.coefficientReal() * 4 == 1,
              "Monomial(\"0.25\") != 0.25");
}
END_TEST

START_TEST (monomial_sum)
{
  mps::formal::Monomial one("1", 0);
  mps::formal::Monomial x("2", 1);

  mps::formal::Polynomial p = one + x;
  
  std::cout << p << std::endl;

  fail_unless (p.degree() == 1, 
	       "deg(2x + 1) != 1");
  fail_unless (p[0].coefficientReal() == 1,
	       "2x+1 has constant term different from 1");
  fail_unless (p[1].coefficientReal() == 2,
	       "2x+1 has leading coefficient different from 2");
}
END_TEST

int
main (void)
{
  int number_failed;

  starting_setup ();

  Suite *s = suite_create ("Formal arithmetic");
  TCase *tc_monomials = tcase_create ("Monomial operations");

  // Basic operations on list elements
  tcase_add_test (tc_monomials, monomial_creation);
  tcase_add_test (tc_monomials, monomial_sum);
  tcase_add_test (tc_monomials, monomial_floating);

  // Basic operations on lists
  suite_add_tcase (s, tc_monomials);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return(number_failed != 0);
}

