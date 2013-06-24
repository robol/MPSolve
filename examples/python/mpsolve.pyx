cdef int n = 4
    # Create a polynomial that will be solved
cdef mps_context *status = mps_context_new ()
cdef mps_monomial_poly *poly = mps_monomial_poly_new (status, n)
# Set the coefficients. We will solve x^n - 1 in here

mps_monomial_poly_set_coefficient_int (status, poly, 0, -1, 0)
mps_monomial_poly_set_coefficient_int (status, poly, n,  1, 0)

# Select some common output options, i.e. 512 bits of precision
# (more or less 200 digits guaranteed) and approximation goal.

mps_context_set_output_prec (status, 512)
mps_context_set_output_goal (status, MPS_OUTPUT_GOAL_APPROXIMATE)

# Solve the polynomial

mps_context_set_input_poly (status, <mps_polynomial*>poly)
mps_mpsolve (status)


# Get the roots in a <code>cplx_t</code> vector. Please note that
# this make completely useless to have asked 512 bits of output
# precision, and you should use mps_context_get_roots_m() to get
# multiprecision approximation of the roots.

cdef cplx_t * results = cplx_valloc (n)
mps_context_get_roots_d (status, &results, NULL)

# Free the data used. This will free the monomial_poly if you have
# not done it by yourself.

for i in range(n):
    print results[i][0]

mps_context_free (status)
cplx_vfree (results)
