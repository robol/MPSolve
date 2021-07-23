#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Try to set some integer/rational coefficients

import mpsolve


def rootsofunity():
    """Solve the polynomial x^n - 1 using rational 
    input for MPSolve"""
    ctx = mpsolve.Context()
    poly = mpsolve.MonomialPoly(ctx, 8)

    poly.set_coefficient (0, "-1")
    poly.set_coefficient (8, "1")

    print( "Roots of x^8 - 1 = %s" % ( "  ".join (map (str, ctx.solve (poly)))))

if __name__ == "__main__":
    rootsofunity()
           

    
