NULL = 

lib_LTLIBRARIES = libmps.la
libmps_la_LDFLAGS = -version-info 3:0:0
includedir=${prefix}/include/mps
include_HEADERS = \
	../../include/mps/tools.h \
	../../include/mps/interface.h \
	../../include/mps/gmptools.h \
	../../include/mps/mpc.h \
	../../include/mps/poly.h \
	../../include/mps/mptemp.h \
	../../include/mps/debug.h \
	../../include/mps/mt.h \
	../../include/mps/link.h \
	../../include/mps/core.h \
	../../include/mps/rursolve.h \
	../../include/mps/secular.h \
	$(NULL)

AM_CPPFLAGS= \
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
	 secular.c \
	 secular-newton.c \
	 secular-starting.c \
	$(NULL)

libmps_la_LIBADD = \
	-lgmp \
	-lm \
	$(NULL)
	
