from mpsolve cimport *

cdef class Context:
    """The Context class is a wrapper around the mps_context type
    in libmps. A Context instance must be instantied before
    allocating and/or solving polynomials and secular equations, 
    and can then be used to specify the desired property of the 
    solution. """
    
    def __init__(self):
        self._c_ctx = mps_context_new()

    def __dealloc__(self):
        if self._c_ctx is not NULL:
            mps_context_free(self._c_ctx)

    def set_input_poly(self, poly):
        self._set_input_poly(poly)

    cdef _set_input_poly(self, Polynomial poly):
        mps_context_set_input_poly(self._c_ctx, <mps_polynomial*> (poly._c_polynomial))

    def mpsolve(self, poly = None):
        if poly is not None:
            self.set_input_poly(poly)
        mps_mpsolve(self._c_ctx)

    def solve(self, poly = None):
        self.mpsolve(poly)
        return self.get_roots()

    def get_roots(self):
        cdef cplx_t * results = NULL
        mps_context_get_roots_d (self._c_ctx, &results, NULL)

        result_list = []
        for i in range(mps_context_get_degree(self._c_ctx)):
            result_list.append(results[i][0])
        return result_list 

    def get_inclusion_radii(self):
        """Return a set of guaranteed inclusion radii for the 
        approximations obtained through a call to get_roots()"""

        cdef cplx_t * results = NULL
        cdef double * radii = NULL
        mps_context_get_roots_d (self._c_ctx, &results, &radii)

        radii_list = []
        for i in range(mps_context_get_degree(self._c_ctx)):
            radii_list.append(radii[i])
        return radii_list
        

cdef class Polynomial:
    """This is a wrapper around mps_polynomial struct. """

    def __dealloc__(self):
        mps_polynomial_free(self._c_ctx, self._c_polynomial)

cdef class MonomialPoly(Polynomial):
    """A polynomial specified with its monomial coefficients. """

    def __cinit__(self, Context ctx, int degree):
        self._degree = degree
        self._c_mp = mps_monomial_poly_new(ctx._c_ctx, degree)
        self._c_polynomial = <mps_polynomial*> self._c_mp
        self._c_ctx = ctx._c_ctx

    def set_coefficient(self, n, coeff):
        """Set coefficient of degree n of the polynomial 
        to the value of coeff. Please note that you should use 
        the same data type for all the coefficients, and you
        should use integers when possible. """

        if n < 0 or n > self._degree:
            raise RuntimeError("Coefficient degree is out of bounds")

        if isinstance(coeff, int):
            mps_monomial_poly_set_coefficient_int(self._c_ctx, 
                                                  self._c_mp, n, coeff, 0)
        elif isinstance(coeff, float):
            mps_monomial_poly_set_coefficient_d(self._c_ctx, self._c_mp, n, coeff, 0)
        else:
            raise RuntimeError("Coefficient type not supported")
