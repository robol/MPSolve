NULL = 

lib_LTLIBRARIES = libmps.la
libmps_la_LDFLAGS = -version-info 3:0:0
includedir=${prefix}/include/mps
include_HEADERS = \
        ${top_srcdir}/include/mps/tools.h \
        ${top_srcdir}/include/mps/interface.h \
	${top_srcdir}/include/mps/input-buffer.h \
        ${top_srcdir}/include/mps/gmptools.h \
        ${top_srcdir}/include/mps/mpc.h \
	${top_srcdir}/include/mps/monomial-poly.h \
	${top_srcdir}/include/mps/cluster.h \
        ${top_srcdir}/include/mps/debug.h \
        ${top_srcdir}/include/mps/mt.h \
        ${top_srcdir}/include/mps/link.h \
        ${top_srcdir}/include/mps/mps.h \
        ${top_srcdir}/include/mps/rursolve.h \
        ${top_srcdir}/include/mps/secular.h \
        ${top_srcdir}/include/mps/threading.h \
	${top_srcdir}/include/mps/status.h \
	${top_srcdir}/include/mps/options.h \
	$(NULL)

AM_CPPFLAGS= \
        -I${top_srcdir}/include \
	$(MPS_CFLAGS) \
	-pthread \
	$(NULL)

libmps_la_SOURCES = \
	    aberth.c \
	    cluster.c \
	    convex.c \
	    data.c \
	    debug.c \
	    defaults.c \
	    general-radius.c \
	    getopts.c \
	    gmptools.c \
	    horner.c \
	    improve.c \
	    inclusion.c \
	    input-buffer.c \
	    interface.c \
	    link.c \
	    main.c \
	    modify.c \
	    monomial-poly.c \
	    monomial-radius.c \
	    monomial-threading.c \
	    mpc.c \
	    mt.c \
	    newton.c \
	    secular-evaluation.c \
	    secular-ga.c \
	    secular-iteration.c \
	    secular-newton.c \
	    secular-radius.c \
	    secular-regeneration.c \
	    secular-starting.c \
	    secular.c \
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

libmps_la_LIBADD = \
	-lgmp \
	-lm \
	$(NULL)
