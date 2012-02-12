NULL = 

bin_PROGRAMS = mpsolve

mpsolve_CFLAGS = \
        -I${top_srcdir}/include \
	$(SECSOLVE_CFLAGS) \
	$(NULL)

mpsolve_SOURCES = \
	mpsolve.c \
	$(NULL)

mpsolve_LDFLAGS = \
	$(SECSOLVE_LIBS)

mpsolve_LDADD = \
	${top_builddir}/src/libmps/libmps.la
