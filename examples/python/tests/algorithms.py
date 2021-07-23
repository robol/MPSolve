#!/usr/bin/env python3
# -*- coding: utf-8
#
# This test shows how you can solve a polynomial using different
# algorithms. 
#
# Author: Leonardo Robol <leonardo.robol@sns.it>
# Date: 13/02/2014

from mpsolve import *

if __name__ == "__main__":

    ctx = Context()
    poly = MonomialPoly(ctx, 4)

    poly.set_coefficient (0, -1)
    poly.set_coefficient (4,  1)

    roots = ctx.solve (poly, Algorithm.SECULAR_GA)
    roots2 = ctx.solve (poly, Algorithm.STANDARD_MPSOLVE)

    print("Polynomial solved by secsolve: %s" % map (complex, roots))
    print("Polynomial solved by unisolve: %s" % map (complex, roots2))
