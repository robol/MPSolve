#!/usr/bin/env python
# -*- coding: utf-8
#
# Example Python script that uses MPSolve to solve various polynomials. 
#
# Author: Leonardo Robol <leonardo.robol@sns.it>
#

import mpsolve

def roots_of_unity(n):
    """Compute the n-th roots of unity, i.e., the solutions
    of the equation x^n = 1"""

    ctx = mpsolve.Context()
    poly = mpsolve.MonomialPoly(ctx, n)

    poly.set_coefficient(n, 1)
    poly.set_coefficient(0, -1)

    return ctx.solve(poly)

if __name__ == "__main__":

    # Get the roots of unity in a fairly simple case
    roots = roots_of_unity (4)
    print("The 4-th roots of unity are: %s" % map (str, roots))

    # Do something more fancy (even if still with the roots
    # of unity). 
    roots = roots_of_unity (1024)

    # You can easily plot the roots using matplotlib, if
    # it's available. We use non interactive plotting here
    # just to make sure that the user can see the plot. 
    from matplotlib import pyplot as plot

    plot.ioff ()
    plot.plot (map (lambda x : x.real, roots),
               map (lambda x : x.imag, roots))
    plot.show ()
    
