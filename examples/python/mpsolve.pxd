cdef extern from "gmp.h":
    cdef struct __mpq_struct:
       pass

    ctypedef __mpq_struct[1] mpq_t

cdef extern from "mps/mps.h":

    cdef struct __mpc_struct:
       pass

    ctypedef __mpc_struct[1] mpc_t

    ctypedef void mps_context
    ctypedef void mps_polynomial
    ctypedef void mps_monomial_poly
    ctypedef void mps_chebyshev_poly
    ctypedef void mps_secular_equation

    ctypedef bint mps_boolean

    cdef enum mps_structure:
      MPS_STRUCTURE_REAL_INTEGER,
      MPS_STRUCTURE_REAL_RATIONAL,
      MPS_STRUCTURE_REAL_FP,
      MPS_STRUCTURE_REAL_BIGFLOAT,
      MPS_STRUCTURE_COMPLEX_INTEGER,
      MPS_STRUCTURE_COMPLEX_RATIONAL,
      MPS_STRUCTURE_COMPLEX_FP,
      MPS_STRUCTURE_COMPLEX_BIGFLOAT,
      MPS_STRUCTURE_UNKNOWN

    mps_polynomial* MPS_POLYNOMIAL (void*)

    # GMP related methods
    cdef void mpq_init(mpq_t)
    cdef void mpq_clear(mpq_t)
    cdef void mpc_init2 (mpc_t, long int)
    cdef void mpc_clear(mpc_t)
    cdef int mpq_set_si(mpq_t, int, int)
    cdef int mpq_set_str(mpq_t, const char *, int)
    cdef void mpc_set_d(mpc_t, double, double)

    # Various helpers
    cdef char* mps_utils_build_equivalent_rational_string(mps_context *, char *)

    #
    # mps_polynomial methods 
    #
    mps_monomial_poly* MPS_MONOMIAL_POLY(void*)
    mps_monomial_poly *mps_monomial_poly_new(mps_context *status, int n)
    void mps_monomial_poly_set_coefficient_int (mps_context * s, mps_monomial_poly * mp, long int i,
                                                long long real_part, long long imag_part)
    void mps_monomial_poly_set_coefficient_d   (mps_context * s, mps_monomial_poly * mp, long int i,
                                                double real_part, double imag_part)
    void mps_monomial_poly_set_coefficient_q   (mps_context * s, mps_monomial_poly * mp, long int i,
                                                mpq_t real_part, mpq_t imag_part)

    #
    # mps_chebyshev_poly methods
    #
    mps_chebyshev_poly* MPS_CHEBYSHEV_POLY(void*)
    mps_chebyshev_poly* mps_chebyshev_poly_new (mps_context * ctx, int n, mps_structure structure)
    void mps_chebyshev_poly_set_coefficient_q (mps_context * ctx, mps_chebyshev_poly *poly, 
                                               int i, mpq_t real_part, mpq_t imag_part)
    void mps_chebyshev_poly_set_coefficient_f (mps_context * ctx, mps_chebyshev_poly * poly, 
                                               int i, mpc_t coeff)

    # 
    # mps_secular_equation methods
    #
    mps_secular_equation* mps_secular_equation_new_raw (mps_context* ctx, int n)    

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

    cdef enum mps_algorithm:
       MPS_ALGORITHM_SECULAR_GA
       MPS_ALGORITHM_STANDARD_MPSOLVE    

    void mps_context_free (mps_context * s)
    void mps_context_select_algorithm (mps_context * ctx, mps_algorithm alg)
    void mps_polynomial_free (mps_context * s, mps_polynomial * poly)


cdef class Context:
    cdef mps_context * _c_ctx
    cdef _set_input_poly(self, Polynomial poly)

cdef class Polynomial:
    cdef mps_polynomial* _c_polynomial
    cdef Context _ctx
    cdef int _degree

cdef class MonomialPoly(Polynomial):
     pass

cdef class ChebyshevPoly(Polynomial):
     pass
