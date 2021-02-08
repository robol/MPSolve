#!/usr/bin/env python3
# -*- coding: utf-8
#
# This tests checks that basic polynomial solving is working from the Python
# module. It is not aimed at detecting all bugs, but should fail in case there
# is something that prevents the Python module to work at all. 
#
# Note that while this catch almost all regressions (the bugs are usually buried
# in the C code that is tested from src/tests) it would be nice to have a complete
# test suite here. 
#
# Author: Leonardo Robol <leonardo.robol@sns.it>
# Date: 13/02/2014

import mpsolve

def simple_polynomial_example():
    """Solve a very simple polynomial. """

    n = 100
    ctx = mpsolve.Context()
    poly = mpsolve.MonomialPoly(ctx, n)

    poly.set_coefficient (0, -1)
    poly.set_coefficient (n,  1)

    print ("x^%d - 1 roots: " % n)
    for root in ctx.solve (poly):
        print (" - " + str(root))

    print("")

def simple_polynomial_int_coefficients():
    n = 2
    ctx = mpsolve.Context()
    poly = mpsolve.MonomialPoly(ctx, n)

    poly.set_coefficient (0, -1)
    poly.set_coefficient (1, 0)
    poly.set_coefficient (2, 2)
    ctx.solve (poly)

    print("")

def simple_polynomial_float_coefficients():
    n = 3
    ctx = mpsolve.Context()
    poly = mpsolve.MonomialPoly(ctx, n)

    poly.set_coefficient (0, -1.0)
    poly.set_coefficient (1, 0.0)
    poly.set_coefficient (2, 2.5)
    ctx.solve (poly)

    print("")

def simple_chebyshev_example():
    """Here is a simple example of a polynomial expressed in the
    Chebyshev basis."""

    n = 100
    ctx = mpsolve.Context()
    poly = mpsolve.ChebyshevPoly(ctx, n)

    poly.set_coefficient (n,  1)

    for root in ctx.solve (poly):
        print (root)
        

if __name__ == "__main__":

    simple_polynomial_example()
    simple_polynomial_int_coefficients()
    simple_polynomial_float_coefficients()
    simple_chebyshev_example()
