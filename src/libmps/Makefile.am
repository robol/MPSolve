NULL = 

lib_LTLIBRARIES = libmps.la
includedir=${prefix}/include/mps
include_HEADERS = \
        ${top_srcdir}/include/mps/approximation.h \
        ${top_srcdir}/include/mps/tools.h \
        ${top_srcdir}/include/mps/interface.h \
        ${top_srcdir}/include/mps/chebyshev.h \
	${top_srcdir}/include/mps/input-buffer.h \
        ${top_srcdir}/include/mps/gmptools.h \
        ${top_srcdir}/include/mps/mpc.h \
	${top_srcdir}/include/mps/monomial-poly.h \
	${top_srcdir}/include/mps/cluster.h \
        ${top_srcdir}/include/mps/debug.h \
	${top_srcdir}/include/mps/mt-types.h \
        ${top_srcdir}/include/mps/mt.h \
        ${top_srcdir}/include/mps/link.h \
        ${top_srcdir}/include/mps/mps.h \
        ${top_srcdir}/include/mps/rursolve.h \
        ${top_srcdir}/include/mps/secular-equation.h \
        ${top_srcdir}/include/mps/threading.h \
	${top_srcdir}/include/mps/context.h \
	${top_srcdir}/include/mps/options.h \
	${top_srcdir}/include/mps/polynomial.h \
	$(NULL)

AM_CPPFLAGS= \
        -I${top_srcdir}/include \
	$(LIBMPS_CFLAGS) \
	$(GMP_CFLAGS) \
	$(PTHREAD_CFLAGS) \
	$(NULL)

libmps_la_SOURCES = \
	    aberth.c \
	    approximation.c \
	    chebyshev.c \
	    chebyshev-evaluation.c \
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
	    horner.c \
	    improve.c \
	    inclusion.c \
	    input-buffer.c \
	    interface.c \
	    jacobi-aberth.c \
	    link.c \
	    main.c \
	    modify.c \
	    monomial-poly.c \
	    monomial-threading.c \
	    mpc.c \
	    mt.c \
	    newton.c \
	    polynomial.c \
	    secular-evaluation.c \
	    secular-ga.c \
	    secular-iteration.c \
	    secular-newton.c \
	    secular-regeneration.c \
	    secular-starting.c \
	    secular-equation.c \
	    solve.c \
	    sort.c \
	    starting.c \
	    stio.c \
	    test.c \
	    threading.c \
	    tools.c \
	    touch.c \
	    user.c \
	    $(NULL)

libmps_la_LDFLAGS = \
        -version-info 3:0:0
	$(PTHREAD_CFLAGS)

libmps_la_LIBADD = \
	$(MPS_LDFLAGS) \
	$(GMP_LIBS) \
	$(PTHREAD_LIBS) \
	$(NULL)
