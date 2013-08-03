cdef extern from "mps/mps.h":

    cdef struct mps_context:
        pass

    cdef struct mps_polynomial:
        pass

    cdef struct mps_monomial_poly:
        pass

    ctypedef bint mps_boolean

    mps_monomial_poly *mps_monomial_poly_new(mps_context *status, int n)
    void mps_monomial_poly_set_coefficient_int (mps_context * s, mps_monomial_poly * mp, long int i,
                                                long long real_part, long long imag_part)
    void mps_monomial_poly_set_coefficient_d   (mps_context * s, mps_monomial_poly * mp, long int i,
                                                double real_part, double imag_part)

    cdef enum mps_output_goal:
        MPS_OUTPUT_GOAL_ISOLATE
        MPS_OUTPUT_GOAL_APPROXIMATE
        MPS_OUTPUT_GOAL_COUNT

    void mps_context_set_output_prec (mps_context * s, long int prec)
    void mps_context_set_input_prec (mps_context * s, long int prec)
    void mps_context_set_output_prec (mps_context * s, long int prec)
    void mps_context_set_output_goal (mps_context * s, mps_output_goal goal)
    void mps_context_set_jacobi_iterations (mps_context * s, mps_boolean jacobi_iterations)
    int mps_context_get_degree (mps_context * s)

    void mps_context_set_input_poly (mps_context * s, mps_polynomial * p)

    ctypedef struct __cplx_struct:
        pass

    ctypedef complex cplx_t[1]

    cplx_t * cplx_valloc(int n)
    int cplx_vfree(cplx_t *V)

    void mps_mpsolve (mps_context * s)

    int mps_context_get_roots_d (mps_context * s, cplx_t ** roots, double **radius)

    mps_context * mps_context_new()

    void mps_context_free (mps_context * s)
    void mps_polynomial_free (mps_context * s, mps_polynomial * poly)


cdef class Context:
    cdef mps_context * _c_ctx
    cdef _set_input_poly(self, Polynomial poly)

cdef class Polynomial:
    cdef mps_polynomial* _c_polynomial
    cdef mps_context* _c_ctx
    cdef int _degree

cdef class MonomialPoly(Polynomial):
    cdef mps_monomial_poly * _c_mp
