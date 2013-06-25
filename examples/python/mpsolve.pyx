cimport mpsolve

cdef class Context:

    cdef mps_context * _c_ctx

    def __init__(self):
        self._c_ctx = mps_context_new()

    def __dealloc__(self):
        if self._c_ctx is not NULL:
            mps_context_free(self._c_ctx)

    def set_input_poly(self, poly):
        poly.load_in_context()

    def mpsolve(self):
        mps_mpsolve(self._c_ctx)
        # mps_output(self._c_ctx)

    def get_roots_d(self):
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

cdef class Polynomial:

    cdef int degree
    cdef mps_context * _c_ctx
    cdef mps_polynomial * cpoly

    def load_in_context(self):
        mps_context_set_input_poly(self._c_ctx, self.cpoly)

cdef class MonomialPoly(Polynomial):

    def __init__(self, degree):
        Polynomial.__init__(self)
        self.degree = degree
        self._c_ctx = NULL
        self.cpoly = NULL

    cdef set_context(self, mps_context * context):
        self._c_ctx = context
        self.cpoly = <mps_polynomial*> mps_monomial_poly_new(context, self.degree)

    def __dealloc__(self):
        pass
        # mps_monomial_poly_free(self._c_ctx, self._c_mp)

    def set_coefficient(self, n, coeff):
        if isinstance(coeff, int):
            mps_monomial_poly_set_coefficient_int(self._c_ctx, 
                                                  <mps_monomial_poly*> self.cpoly, 
                                                  n, coeff, 0)
        else:
            raise RuntimeError("Coefficient type not supported")
