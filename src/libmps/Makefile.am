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

BUILT_SOURCES = monomial/yacc-parser.h
AM_YFLAGS = -d

libmps_la_SOURCES = \
	chebyshev/chebyshev-evaluation.c \
	chebyshev/chebyshev-parser.c \
	chebyshev/chebyshev.c \
	common/aberth.c \
	common/approximation.c \
	common/cluster-analysis.c \
	common/cluster.c \
	common/context.c \
	common/convex.c \
	common/defaults.c \
	common/file-starting.c \
	common/improve.c \
	common/inclusion.c \
	common/inline-poly-parser.c \
	common/interface.c \
	common/jacobi-aberth.c \
	common/list-element.c \
	common/list.c \
	common/modify.c \
	common/nroots-polynomial.cpp \
	common/parser.c \
	common/polynomial.c \
	common/polynomialxx.cpp \
	common/recursive-starting.c \
	common/sort.c \
	common/starting-configuration.c \
	common/starting.c \
        common/strndup.c \
	common/test.c \
	common/tools.c \
	common/touch.c \
	common/user.c \
	common/utils.c \
	common/version.c \
	common/validation.c \
	formal/formal-monomial.cpp \
	formal/formal-polynomial.cpp \
	floating-point/gmptools.c \
	floating-point/link.c \
	floating-point/mpcf.c \
	floating-point/mt.c \
	general/general-radius.c \
	general/general-starting.c \
	matrix/hessenberg-determinant.c \
	monomial/horner.c \
	monomial/monomial-matrix-poly.c \
	monomial/monomial-parser.c \
	monomial/monomial-poly.c \
	monomial/monomial-threading.c \
	monomial/newton.c \
	monomial/yacc-parser.y \
	monomial/tokenizer.l \
	monomial/shift.c \
	secsolve/secular-ga.c \
	secsolve/secular-iteration.c \
	secsolve/secular-regeneration.c \
	secsolve/standard-regeneration-driver.c \
	secular/secular-equation.c \
	secular/secular-evaluation.c \
	secular/secular-newton.c \
	secular/secular-parser.c \
	secular/secular-starting.c \
	system/abstract-input-stream.cpp \
	system/file-input-stream.cpp \
	system/memory-file-stream.cpp \
	system/data.c \
	system/debug.c \
	system/getline.c \
	system/getopts.c \
	system/input-buffer.c \
	system/input-output.c \
	system/threading.c \
	unisolve/main.c \
	unisolve/solve.c \
	$(NULL)

libmps_la_LDFLAGS = \
        -version-info 3:1:0 \
	-no-undefined \
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
	$(PTHREAD_CFLAGS) \
	-no-undefined

libmpsprivate_la_SOURCES = \
	$(libmps_la_SOURCES)

libmpsprivate_la_LIBADD = \
	$(libmps_la_LIBADD)

libmpsprivate_la_CFLAGS = \
	-DMPS_PUBLISH_PRIVATE_METHODS=1

endif
