import mpsolve

n = 1200

ctx = mpsolve.Context()
poly = ctx.monomial_poly_create(n)

poly.set_coefficient(0, -1)
poly.set_coefficient(n, 1)

ctx.set_input_poly(poly)
ctx.mpsolve()

roots = ctx.get_roots_d()
print "\n".join(map(str, roots))
