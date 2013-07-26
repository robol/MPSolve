from mpsolve cimport *

cdef class Context:

    def __init__(self):
        self._c_ctx = mps_context_new()

    def __dealloc__(self):
        if self._c_ctx is not NULL:
            mps_context_free(self._c_ctx)

    def set_input_poly(self, poly):
        self._set_input_poly(poly)

    cdef _set_input_poly(self, Polynomial poly):
        mps_context_set_input_poly(self._c_ctx, <mps_polynomial*> (poly._c_polynomial))

    def mpsolve(self):
        mps_mpsolve(self._c_ctx)

    def get_roots(self):
        cdef cplx_t * results = NULL
        mps_context_get_roots_d (self._c_ctx, &results, NULL)

        result_list = []
        for i in range(mps_context_get_degree(self._c_ctx)):
            result_list.append(results[i][0])
        return result_list        

cdef class Polynomial:
    def __dealloc__(self):
        mps_polynomial_free(self._c_ctx, self._c_polynomial)



cdef class MonomialPoly:

    def __cinit__(self, Context ctx, int degree):
        self._degree = degree
        self._c_mp = mps_monomial_poly_new(ctx._c_ctx, degree)
        self._c_polynomial = <mps_polynomial*> self._c_mp
        self._c_ctx = ctx._c_ctx

    def set_coefficient(self, n, coeff):
        if isinstance(coeff, int):
            mps_monomial_poly_set_coefficient_int(self._c_ctx, 
                                                  self._c_mp, n, coeff, 0)
        else:
            raise RuntimeError("Coefficient type not supported")
