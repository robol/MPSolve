NULL = 

bin_PROGRAMS = unisolve

unisolve_CFLAGS = \
        -I${top_srcdir}/include \
	$(UNISOLVE_CFLAGS) \
	$(NULL)

unisolve_SOURCES = \
	unisolve.c \
	$(NULL)

unisolve_LDFLAGS = \
	$(UNISOLVE_LIBS)
	$(NULL)

unisolve_LDADD = \
	$(top_builddir)/src/libmps/libmps.la \
	$(NULL)

unisolve_test_LDADD = \
	$(top_builddir)/src/libmps/libmps.la \
	$(NULL)
