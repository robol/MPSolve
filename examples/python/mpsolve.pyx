from mpsolve cimport *

cdef class Context:

    def __init__(self):
        self._c_ctx = mps_context_new()

    def __dealloc__(self):
        if self._c_ctx is not NULL:
            mps_context_free(self._c_ctx)

    def set_input_poly(self, poly):
        self._set_input_poly(poly)

    cdef _set_input_poly(self, MonomialPoly poly):
        mps_context_set_input_poly(self._c_ctx, <mps_polynomial*> (poly._c_mp))
        # poly.load_into_context(self._c_ctx)

    def mpsolve(self):
        mps_mpsolve(self._c_ctx)

    def get_roots(self):
        cdef cplx_t * results = NULL
        mps_context_get_roots_d (self._c_ctx, &results, NULL)

        result_list = []
        for i in range(mps_context_get_degree(self._c_ctx)):
            result_list.append(results[i][0])
        return result_list

    def monomial_poly_create(self, degree):
        poly = MonomialPoly(degree)
        poly.set_context(self._c_ctx)
        return poly
        

cdef class MonomialPoly:

    # cdef mps_monomial_poly * _c_mp

    def __init__(self, degree):
        self._degree = degree
        self._c_mp = NULL
        self._c_ctx = NULL

    cdef load_into_context(self, mps_context * context):
        mps_context_set_input_poly(context, <mps_polynomial*> self._c_mp)

    cdef set_context(self, mps_context * context):
        self._c_ctx = context
        self._c_mp = mps_monomial_poly_new(context, self._degree)

    def __dealloc__(self):
        pass
        # mps_monomial_poly_free(self._c_ctx, self._c_mp)

    def set_coefficient(self, n, coeff):
        if isinstance(coeff, int):
            mps_monomial_poly_set_coefficient_int(self._c_ctx, 
                                                  self._c_mp, n, coeff, 0)
        else:
            raise RuntimeError("Coefficient type not supported")
