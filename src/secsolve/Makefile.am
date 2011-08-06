NULL = 

bin_PROGRAMS = secsolve

secsolve_CFLAGS = \
        -I${top_srcdir}/include \
	$(SECSOLVE_CFLAGS) \
	$(NULL)

secsolve_SOURCES = \
	secsolve.c \
	$(NULL)

secsolve_LDFLAGS = \
	$(SECSOLVE_LIBS)

secsolve_LDADD = \
	${top_builddir}/src/libmps/libmps.la
