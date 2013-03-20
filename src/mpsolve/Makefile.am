NULL = 

bin_PROGRAMS = mpsolve

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
