NULL = 

bin_PROGRAMS = mpsolve

man_MANS = ${srcdir}/mpsolve.1

EXTRA_DIST = mpsolve.1

mpsolve_CFLAGS = \
        -I${top_srcdir}/include \
	-I. \
	$(MPSOLVE_CFLAGS) \
	$(GTK_CFLAGS) \
	$(NULL)

mpsolve_SOURCES = \
	iteration-logger.h \
	iteration-logger.c \
	mpsolve.c \
	$(NULL)

mpsolve_LDFLAGS = \
	$(MPSOLVE_LIBS) \
	$(GTK_LIBS)

mpsolve_LDADD = \
	${top_builddir}/src/libmps/libmps.la
