from mpsolve cimport *

class Algorithm:
    """Here you can find all the available algorithms in MPSolve. 
    You can use these contants to specify the algorithm of your choice
    when you call Context.solve() or Context.mpsolve(). For example: 
    
     poly = MonomialPoly(ctx, n)
     roots = ctx.solve (poly, Algorithm.SECULAR_GA) 
    """
  
    SECULAR_GA = MPS_ALGORITHM_SECULAR_GA
    STANDARD_MPSOLVE = MPS_ALGORITHM_STANDARD_MPSOLVE

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
        """Select the polynomial that should be solved when mpsolve() is called.
        Note that each Context can only solve one polynomial at a time."""
        self._set_input_poly(poly)

    cdef _set_input_poly(self, Polynomial poly):
        mps_context_set_input_poly(self._c_ctx, <mps_polynomial*> (poly._c_polynomial))

    def mpsolve(self, poly = None, algorithm = MPS_ALGORITHM_SECULAR_GA):
        """Calling this method will trigger the solution of the polynomial
        previously loaded by a call to set_input_poly, or to the one passed
        as second argument to this function. 
        
        An optional third argument specify the desired algorithm. """
        if poly is not None:
            self.set_input_poly(poly)

        # Select the proper algorithm for this polynomial
        if isinstance (poly, ChebyshevPoly):
            algorithm = MPS_ALGORITHM_SECULAR_GA

        mps_context_select_algorithm (self._c_ctx, algorithm)

        mps_mpsolve(self._c_ctx)

    def solve(self, poly = None, algorithm = MPS_ALGORITHM_SECULAR_GA):
        """Simple shorthand for the combination of set_input_poly() and mpsolve(). 
        This function directly returns the approximations that could otherwise be
        obtained by a call to the get_roots() method. """
        self.mpsolve(poly)
        return self.get_roots()

    def get_roots(self):
        """Returns  the approximations obtained by MPSolve after a call to the mpsolve
        method. Consider using the convienience solve() method, that is usually more
        convenient."""
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

    def __cinit__(self, Context ctx, int degree):
        self._degree = degree
        self._ctx = ctx

    def __dealloc__(self):
        mps_polynomial_free(self._ctx._c_ctx, self._c_polynomial)

cdef class MonomialPoly(Polynomial):
    """A polynomial specified with its monomial coefficients. """

    def __cinit__(self, Context ctx, int degree):
        Polynomial.__init__(self, ctx, degree)
        self._c_polynomial = MPS_POLYNOMIAL (mps_monomial_poly_new (ctx._c_ctx, degree))

    def set_coefficient(self, n, coeff):
        """Set coefficient of degree n of the polynomial 
        to the value of coeff. Please note that you should use 
        the same data type for all the coefficients, and you
        should use integers when possible. """

        cdef mps_monomial_poly *mp = MPS_MONOMIAL_POLY (self._c_polynomial)
        cdef mpq_t qcoeff, qzero

        mpq_init (qcoeff)
        mpq_init (qzero)

        if n < 0 or n > self._degree:
            raise RuntimeError("Coefficient degree is out of bounds")

        if isinstance(coeff, int):
            mps_monomial_poly_set_coefficient_int(self._ctx._c_ctx, 
                                                  mp, n, coeff, 0)
        elif isinstance(coeff, float):
            mps_monomial_poly_set_coefficient_d(self._ctx._c_ctx, mp, n, coeff, 0)
        elif isinstance(coeff, str):
            # We expect the user to specify the input as string in case of rational
            # input. 
            eq = mps_utils_build_equivalent_rational_string (self._ctx._c_ctx, 
                                                             coeff)
            mpq_set_str (qcoeff, eq, 10)
            mps_monomial_poly_set_coefficient_q (self._ctx._c_ctx, 
                                                 mp, n, qcoeff, qzero)

        else:
            mpq_clear (qzero)
            mpq_clear (qcoeff)
            raise RuntimeError("Coefficient type not supported")

        mpq_clear (qzero)
        mpq_clear (qcoeff)


cdef class ChebyshevPoly(Polynomial):
    """A polynomial represented in the Chebyshev base."""

    def __cinit__(self, Context ctx, int degree):
        Polynomial.__init__(self, ctx, degree)
        self._c_polynomial = MPS_POLYNOMIAL (
            mps_chebyshev_poly_new (ctx._c_ctx, 
                                    degree, 
                                    MPS_STRUCTURE_COMPLEX_RATIONAL))

    def set_coefficient(self, n, coeff):
        """Set the coefficient of degree n of the polynomial"""
        cdef mps_chebyshev_poly* cb = MPS_CHEBYSHEV_POLY (self._c_polynomial)

        if n < 0 or n > self._degree:
            raise RuntimeError("Coefficient degree is out of bounds")

        cdef mpq_t qcoeff
        cdef mpq_t qzero
        cdef mpc_t ccoeff

        cdef char * eq

        mpq_init (qcoeff)
        mpq_init (qzero)
        mpc_init2 (ccoeff, 64)

        mpq_set_si (qzero, 0, 1)

        if isinstance(coeff, int):
            mpq_set_si (qcoeff, coeff, 1)
            mps_chebyshev_poly_set_coefficient_q (self._ctx._c_ctx, 
                                                  cb, n, qcoeff, qzero)
        elif isinstance(coeff, float):
            mpc_set_d (ccoeff, coeff, 0.0)
        elif isinstance (coeff, str):
            eq = mps_utils_build_equivalent_rational_string (self._ctx._c_ctx, 
                                                             coeff)
            mpq_set_str (qcoeff, eq, 10)
            mps_chebyshev_poly_set_coefficient_q (self._ctx._c_ctx, 
                                                  cb, n, qcoeff, qzero)
        else:
            mpq_clear (qcoeff)
            mpq_clear (qzero)
            mpc_clear (ccoeff)
            raise RuntimeError("Unsupported type for coefficient")

        mpq_clear (qcoeff)
        mpc_clear (ccoeff)
        mpq_clear (qzero)

        
            
