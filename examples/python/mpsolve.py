import sys
import ctypes
import ctypes.util


# Load the libmps shared library. We should keep the .so version update
# in case we bump it in the future.
_mps = ctypes.CDLL("libmps.so.3")


class Cplx(ctypes.Structure):
    """Wrapper around the cplx_t type of MPSolve, that is usually
    a direct mapping onto the complex type of C99, but has a fallback
    custom implementation for systems that do not provide the type"""

    _fields_ = [("real", ctypes.c_double),
                ("imag", ctypes.c_double)]

    def __repr__(self):
        return "%e + %ei" % (self.real, self.imag)

    def __complex__(self):
        return complex(self.real + 1j*self.imag)


class Goal:
    """ Goal to reach before returning the result.
    """
    MPS_OUTPUT_GOAL_ISOLATE = 0
    MPS_OUTPUT_GOAL_APPROXIMATE = 1
    MPS_OUTPUT_GOAL_COUNT = 2


class Algorithm:
    """Here you can find all the available algorithms in MPSolve.
    You can use these contants to specify the algorithm of your choice
    when you call Context.solve() or Context.mpsolve(). For example:

     poly = MonomialPoly(ctx, n)
     roots = ctx.solve(poly, Algorithm.SECULAR_GA)
    """

    STANDARD_MPSOLVE = 0
    SECULAR_GA = 1


class Context:
    """The Context class is a wrapper around the mps_context type
    in libmps. A Context instance must be instantied before
    allocating and/or solving polynomials and secular equations,
    and can then be used to specify the desired property of the
    solution. """

    def __init__(self):
        self._c_ctx = _mps.mps_context_new()

    def __del__(self):
        if self._c_ctx is not None:
            _mps.mps_context_free(self._c_ctx)

    def set_input_poly(self, poly):
        """Select the polynomial that should be solved when mpsolve() is
        called. Note that each Context can only solve one polynomial
        at a time."""
        self._set_input_poly(poly)

    def _set_input_poly(self, poly):
        _mps.mps_context_set_input_poly(self._c_ctx, poly._c_polynomial)

    def mpsolve(self, poly=None, algorithm=Algorithm.SECULAR_GA):
        """Calling this method will trigger the solution of the polynomial
        previously loaded by a call to set_input_poly, or to the one passed
        as second argument to this function.

        An optional third argument specify the desired algorithm. """
        if poly is not None:
            self.set_input_poly(poly)

        # Select the proper algorithm for this polynomial
        if isinstance(poly, ChebyshevPoly):
            algorithm = Algorithm.SECULAR_GA

        _mps.mps_context_select_algorithm(self._c_ctx, algorithm)
        # _mps.mps_context_set_output_prec(self._c_ctx,
        #                                 Goal.MPS_OUTPUT_GOAL_APPROXIMATE)
        # _mps.mps_context_set_input_prec(self._c_ctx, ctypes.c_long(0))
        # _mps.mps_context_set_output_prec(self._c_ctx, ctypes.c_long(60))
        _mps.mps_mpsolve(self._c_ctx)

    def solve(self, poly=None, algorithm=Algorithm.SECULAR_GA):
        """Simple shorthand for the combination of set_input_poly() and
        mpsolve(). This function directly returns the approximations
        that could otherwise be obtained by a call to the get_roots() method.
        """
        self.mpsolve(poly)
        return self.get_roots()

    def get_roots(self):
        """Returns  the approximations obtained by MPSolve after a call
        to the mpsolve method. Consider using the convienience solve() method
        instead."""
        degree = _mps.mps_context_get_degree(self._c_ctx)
        results = (Cplx*degree)()

        _mps.mps_context_get_roots_d(self._c_ctx,
                                     ctypes.pointer(ctypes.pointer(results)),
                                     None)
        return [complex(x) for x in results]

    def get_inclusion_radii(self):
        """Return a set of guaranteed inclusion radii for the
        approximations obtained through a call to get_roots()"""
        degree = _mps.mps_context_get_degree(self._c_ctx)
        results = (Cplx*degree)()
        radii = (ctypes.c_double*degree)()

        _mps.mps_context_get_roots_d(self._c_ctx,
                                     ctypes.pointer(ctypes.pointer(results)),
                                     ctypes.pointer(ctypes.pointer(radii)))

        return [float(x) for x in radii]


class Polynomial:
    """This is a wrapper around mps_polynomial struct. """

    def __init__(self, ctx, degree):
        self._degree = int(degree)
        self._ctx = ctx

    def __del__(self):
        _mps.mps_polynomial_free(self._ctx._c_ctx, self._c_polynomial)


class MonomialPoly(Polynomial):
    """A polynomial specified with its monomial coefficients. """

    def __init__(self, ctx, degree):
        Polynomial.__init__(self, ctx, degree)
        deg_c = ctypes.c_long(degree)
        self._c_polynomial = _mps.mps_monomial_poly_new(ctx._c_ctx, deg_c)

    def set_coefficient(self, n, coeff_re, coeff_im=None):
        """Set coefficient of degree n of the polynomial
        to the value of coeff. Please note that you should use
        the same data type for all the coefficients, and you
        should use integers when possible. """

        if coeff_im is not None and type(coeff_re) != type(coeff_im):
            raise ValueError("Coefficient's real and imaginary parts \
have different types")

        mp = self._c_polynomial
        cntxt = self._ctx._c_ctx
        if n < 0 or n > self._degree:
            raise RuntimeError("Coefficient degree is out of bounds")

        if isinstance(coeff_re, int):
            if coeff_im is None:
                coeff_im = 0
            coeff_re = ctypes.c_longlong(coeff_re)
            coeff_im = ctypes.c_longlong(coeff_im)
            _mps.mps_monomial_poly_set_coefficient_int(cntxt, mp, n,
                                                       coeff_re, coeff_im)
        elif isinstance(coeff_re, float):
            if coeff_im is None:
                coeff_im = 0.0
            coeff_re = ctypes.c_double(coeff_re)
            coeff_im = ctypes.c_double(coeff_im)
            _mps.mps_monomial_poly_set_coefficient_d(cntxt, mp, n,
                                                     coeff_re, coeff_im)
        elif isinstance(coeff_re, str):
            if coeff_im is None:
                coeff_im = "0.0"
            if sys.version_info.major > 2:
                coeff_re = bytes(coeff_re, "ASCII")
                coeff_im = bytes(coeff_im, "ASCII")
            _mps.mps_monomial_poly_set_coefficient_s(cntxt, mp, n,
                                                     coeff_re, coeff_im)
        else:
            raise RuntimeError("Coefficient type not supported")

    def get_coefficient(self, n):
        """ Get a coefficient of the polynomial
        """
        mp = self._c_polynomial
        if n < 0 or n > self._degree or not isinstance(n, int):
            raise ValueError("Invalid coefficient degree")

        cf = (Cplx)()
        _mps.mps_monomial_poly_get_coefficient_d(self._ctx._c_ctx, mp, n,
                                                 ctypes.pointer(cf))
        return complex(cf)

    def get_coefficients(self):
        """ Get list of coefficients of the polynomial
        """
        mp = self._c_polynomial
        cf = (Cplx)()
        coeffs = []
        for n in range(self._degree + 1):
            _mps.mps_monomial_poly_get_coefficient_d(self._ctx._c_ctx, mp, n,
                                                     ctypes.pointer(cf))
            cf_n = complex(cf)
            coeffs.append(cf_n)
        return coeffs


class ChebyshevPoly(Polynomial):
    """A polynomial represented in the Chebyshev base."""

    def __init__(self, ctx, degree):
        Polynomial.__init__(self, ctx, degree)

        # 5 here is the equivalent of MPS_STRUCTURE_COMPLEX_RATIONAL
        self._c_polynomial = _mps.mps_chebyshev_poly_new(ctx._c_ctx,
                                                         int(degree), 5)

    def set_coefficient(self, n, coeff_re, coeff_im=None):
        """Set the coefficient of degree n of the polynomial"""
        cb = self._c_polynomial

        if coeff_im is not None and type(coeff_re) != type(coeff_im):
            raise ValueError("Coefficient's real and imaginary parts \
have different types")

        mp = self._c_polynomial
        cntxt = self._ctx._c_ctx

        if n < 0 or n > self._degree:
            raise RuntimeError("Coefficient degree is out of bounds")

        if isinstance(coeff_re, int):
            if coeff_im is None:
                coeff_im = 0
            _mps.mps_chebyshev_poly_set_coefficient_i(cntxt, cb, n,
                                                      coeff_re, coeff_im)
        elif isinstance(coeff, float):
            if coeff_im is None:
                coeff_im = 0.0
            _mps.mpc_set_d(ccoeff, coeff_re, coeff_im)
        # elif isinstance(coeff, str):
        #     if coeff_im is None:
        #         coeff_im = "0.0"
        #     _mps.mps_chebyshev_poly_set_coefficient_s(cntxt, cb, n,
        #                                               coeff_re, coeff_im)
        else:
            raise RuntimeError("Unsupported type for coefficient")
