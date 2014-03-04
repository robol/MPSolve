NULL = 

lib_LTLIBRARIES = libmps.la

includedir=${prefix}/include/mps

AM_CPPFLAGS= \
        -I${top_srcdir}/include \
        -I${top_builddir}/include \
	$(LIBMPS_CFLAGS) \
	$(GMP_CFLAGS) \
	$(PTHREAD_CFLAGS) \
	$(NULL)

libmps_la_SOURCES = \
	    aberth.c \
	    approximation.c \
	    chebyshev.c \
	    chebyshev-evaluation.c \
	    chebyshev-parser.c \
	    cluster.c \
	    cluster-analysis.c \
	    context.c \
	    convex.c \
	    data.c \
	    debug.c \
	    defaults.c \
	    general-radius.c \
	    general-starting.c \
	    getopts.c \
	    gmptools.c \
	    hessenberg-determinant.c \
	    horner.c \
	    improve.c \
	    inclusion.c \
	    inline-poly-parser.c \
	    input-buffer.c \
	    input-output.c \
	    interface.c \
	    jacobi-aberth.c \
	    link.c \
	    main.c \
	    modify.c \
	    monomial-parser.c \
	    monomial-matrix-poly.c \
	    monomial-poly.c \
	    monomial-threading.c \
	    mpc.c \
	    mt.c \
	    newton.c \
	    parser.c \
	    polynomial.c \
	    recursive-starting.c \
	    secular-equation.c \
	    secular-evaluation.c \
	    secular-ga.c \
	    secular-iteration.c \
	    secular-parser.c \
	    secular-newton.c \
	    secular-regeneration.c \
	    secular-starting.c \
	    shift.c \
	    solve.c \
	    sort.c \
	    starting.c \
	    test.c \
	    threading.c \
	    tools.c \
	    touch.c \
	    user.c \
	    utils.c \
	    validation.c \
	    $(NULL)

libmps_la_LDFLAGS = \
        -version-info 3:0:0 \
	$(PTHREAD_CFLAGS)

libmps_la_LIBADD = \
	$(MPS_LDFLAGS) \
	$(GMP_LIBS) \
	$(PTHREAD_LIBS) \
	$(NULL)

#
# In case we need to run the test-suite compile a second version
# of the mps library with all the routines made public, and so, testable :)
#
if CHECK

check_LTLIBRARIES = libmpsprivate.la 	

libmpsprivate_la_CPPFLAGS = \
        -I${top_srcdir}/include \
        -I${top_builddir}/include \
	$(LIBMPS_CFLAGS) \
	$(GMP_CFLAGS) \
	$(PTHREAD_CFLAGS) \
	$(NULL)

libmpsprivate_la_LDFLAGS = \
	$(PTHREAD_CFLAGS)

libmpsprivate_la_SOURCES = \
	$(libmps_la_SOURCES)

libmpsprivate_la_LIBADD = \
	$(libmps_la_LIBADD)

libmpsprivate_la_CFLAGS = \
	-DMPS_PUBLISH_PRIVATE_METHODS=1

endif
