import mpsolve as mps
from time import time
from sys import argv
#
if 'M' not in locals():
    if len(argv) < 2:
        M = 20
    else:
        M = int(argv[1])

phi = [1]
for m in range(M):
    phi0 = list(phi)
    phi = [0] + phi0
    phi[:m+1] = [x - (m + 1)*y for x, y in zip(phi, phi0)]
del phi0
ctx = mps.Context()
poly = mps.MonomialPoly(ctx, M)
for k in range(len(phi)):
    poly.set_coefficient(k, "{:d}".format(phi[k]))
cf = poly.get_coefficients()
# cf_err = [abs(x - y) for x,y in zip(cf, phi)]
# print("\nget_coefficient() error: {}\n".format(max(cf_err)))
t_begin = time()
roots = ctx.solve(poly)
t_end = time()
roots = sorted(roots, key=lambda x: x.real)
err = [abs(x - y - 1) for x, y in zip(roots, range(M))]
radii = ctx.get_inclusion_radii()
print("Roots:\t\t\t\t\tError\t\tInclusion")
print("real part\t imaginary part\t\t\t\tradius")
print("="*65)
for m in range(M):
    print("{:.3e}\t{:.3e}\t\t{:.3e}\t{:.3e}"
          .format(roots[m].real, roots[m].imag, err[m], radii[m]))
print("="*65)
print("\t\t\t\tmax(Error) = {:.3e}".format(max(err)))
t_elapsed = t_end - t_begin
print("\n\tRoot finding took {:d} min {:.3f} sec".
      format(int(t_elapsed//60), t_elapsed % 60))
