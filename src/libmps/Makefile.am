NULL = 

lib_LTLIBRARIES = libmps.la
libmps_la_LDFLAGS = -version-info 3:0:0
includedir=${prefix}/include/mps
include_HEADERS = \
        ${top_srcdir}/include/mps/tools.h \
        ${top_srcdir}/include/mps/interface.h \
        ${top_srcdir}/include/mps/gmptools.h \
        ${top_srcdir}/include/mps/mpc.h \
        ${top_srcdir}/include/mps/poly.h \
        ${top_srcdir}/include/mps/mptemp.h \
        ${top_srcdir}/include/mps/debug.h \
        ${top_srcdir}/include/mps/mt.h \
        ${top_srcdir}/include/mps/link.h \
        ${top_srcdir}/include/mps/core.h \
        ${top_srcdir}/include/mps/rursolve.h \
        ${top_srcdir}/include/mps/secular.h \
        ${top_srcdir}/include/mps/threading.h \
	$(NULL)

AM_CPPFLAGS= \
        -I${top_srcdir}/include \
	$(MPS_CFLAGS) \
	-pthread \
	$(NULL)

libmps_la_SOURCES = \
	 user.c \
	 getopts.c \
	 solve.c \
	 main.c \
	 data.c \
	 cluster.c \
	 mptemp.c \
	 stio.c \
	 gmptools.c \
	 defaults.c \
	 sort.c \
	 tools.c \
	 poly.c \
	 starting.c \
	 convex.c \
	 mt.c \
	 aberth.c \
	 improve.c \
	 test.c \
	 newton.c \
	 mpc.c \
	 link.c \
	 touch.c \
	 interface.c \
	 threading.c \
	 secular.c \
	 secular-newton.c \
	 secular-starting.c \
	 secular-ga.c \
	$(NULL)

libmps_la_LIBADD = \
	-lgmp \
	-lm \
	$(NULL)
	
