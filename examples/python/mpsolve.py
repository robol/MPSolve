import ctypes
import ctypes.util

# Load the libmps shared library. We should keep the .so version update
# in case we bump it in the future. 
_mps = ctypes.CDLL ("libmps.so.3")

_mps.mps_chebyshev_poly_new.restype = ctypes.c_void_p
_mps.mps_context_new.restype = ctypes.c_void_p
_mps.mps_monomial_poly_new.restype = ctypes.c_void_p

class Cplx (ctypes.Structure):
    """Wrapper around the cplx_t type of MPSolve, that is usually
    a direct mapping onto the complex type of C99, but has a fallback
    custom implementation for systems that do not provide the type"""

    _fields_ = [ ("real", ctypes.c_double),
                 ("imag", ctypes.c_double) ]

    def __repr__(self):
        return "%e + %ei" % (self.real, self.imag)

    def __complex__(self):
        return complex(self.real + 1j * self.imag)

class Algorithm:
    """Here you can find all the available algorithms in MPSolve. 
    You can use these contants to specify the algorithm of your choice
    when you call Context.solve() or Context.mpsolve(). For example: 
    
     poly = MonomialPoly(ctx, n)
     roots = ctx.solve (poly, Algorithm.SECULAR_GA) 
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
        self._c_ctx = ctypes.c_void_p(_mps.mps_context_new())

    def __del__(self):
        if self._c_ctx is not None:
            _mps.mps_context_free(self._c_ctx)

    def set_input_poly(self, poly):
        """Select the polynomial that should be solved when mpsolve() is called.
        Note that each Context can only solve one polynomial at a time."""
        self._set_input_poly(poly)

    def _set_input_poly(self, poly):
        _mps.mps_context_set_input_poly(self._c_ctx, poly._c_polynomial)

    def mpsolve(self, poly = None, algorithm = Algorithm.SECULAR_GA):
        """Calling this method will trigger the solution of the polynomial
        previously loaded by a call to set_input_poly, or to the one passed
        as second argument to this function. 
        
        An optional third argument specify the desired algorithm. """
        if poly is not None:
            self.set_input_poly(poly)

        # Select the proper algorithm for this polynomial
        if isinstance (poly, ChebyshevPoly):
            algorithm = Algorithm.SECULAR_GA

        _mps.mps_context_select_algorithm (self._c_ctx, algorithm)
        _mps.mps_mpsolve(self._c_ctx)

    def solve(self, poly = None, algorithm = Algorithm.SECULAR_GA):
        """Simple shorthand for the combination of set_input_poly() and mpsolve(). 
        This function directly returns the approximations that could otherwise be
        obtained by a call to the get_roots() method. """
        self.mpsolve(poly, algorithm)
        return self.get_roots()

    def get_roots(self):
        """Returns  the approximations obtained by MPSolve after a call to the mpsolve
        method. Consider using the convienience solve() method, that is usually more
        convenient."""
        degree = _mps.mps_context_get_degree (self._c_ctx)
        results = (Cplx * degree)()

        _mps.mps_context_get_roots_d (self._c_ctx, 
                                        ctypes.pointer(ctypes.pointer(results)), 
                                        None)

        return results

    def get_inclusion_radii(self):
        """Return a set of guaranteed inclusion radii for the 
        approximations obtained through a call to get_roots()"""
        degree = _mps.mps_context_get_degree (self._c_ctx)
        results = (Cplx * degree)()
        radii = (ctypes.c_double * degree)()

        _mps.mps_context_get_roots_d (self._c_ctx,
                                        ctypes.pointer(ctypes.pointer(results)),
                                        ctypes.pointer(ctypes.pointer(radii)))

        return radii

class Polynomial:
    """This is a wrapper around mps_polynomial struct. """

    def __init__(self, ctx, degree):
        self._degree = degree
        self._ctx = ctx

    def __del__(self):
        _mps.mps_polynomial_free(self._ctx._c_ctx, self._c_polynomial)

class MonomialPoly(Polynomial):
    """A polynomial specified with its monomial coefficients. """

    def __init__(self, ctx, degree):
        Polynomial.__init__(self, ctx, degree)
        self._c_polynomial = \
            ctypes.c_void_p(_mps.mps_monomial_poly_new (ctx._c_ctx, degree))

    def set_coefficient(self, n, coeff):
        """Set coefficient of degree n of the polynomial 
        to the value of coeff. Please note that you should use 
        the same data type for all the coefficients, and you
        should use integers when possible. """

        mp = self._c_polynomial

        if n < 0 or n > self._degree:
            raise RuntimeError("Coefficient degree is out of bounds")

        if isinstance(coeff, int):
            _mps.mps_monomial_poly_set_coefficient_int(self._ctx._c_ctx, 
                                                         mp, n, coeff, 0)
        elif isinstance(coeff, float):
            _mps.mps_monomial_poly_set_coefficient_d(self._ctx._c_ctx, mp, n, coeff, 0)
        elif isinstance(coeff, str):
            _mps.mps_monomial_poly_set_coefficient_s(self._ctx._c_ctx, mp, n, coeff, None);
        else:
            raise RuntimeError("Coefficient type not supported")


class ChebyshevPoly(Polynomial):
    """A polynomial represented in the Chebyshev base."""

    def __init__(self, ctx, degree):
        Polynomial.__init__(self, ctx, degree)

        # 5 here is the equivalent of MPS_STRUCTURE_COMPLEX_RATIONAL
        self._c_polynomial = \
            ctypes.c_void_p(_mps.mps_chebyshev_poly_new (ctx._c_ctx,
                                                          degree, 5))

    def set_coefficient(self, n, coeff):
        """Set the coefficient of degree n of the polynomial"""
        cb = self._c_polynomial

        if n < 0 or n > self._degree:
            raise RuntimeError("Coefficient degree is out of bounds")

        if isinstance(coeff, int):
            _mps.mps_chebyshev_poly_set_coefficient_i (self._ctx._c_ctx, 
                                                       cb, n, coeff, 0)
        elif isinstance(coeff, float):
            _mps.mpc_set_d (ccoeff, coeff, 0.0)
        elif isinstance (coeff, str):
            _mps.mps_chebyshev_poly_set_coefficient_s (self._ctx._c_ctx, 
                                                       cb, n, coeff, "0.0")
        else:
            raise RuntimeError("Unsupported type for coefficient")
        
            
