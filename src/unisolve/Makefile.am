NULL = 

bin_PROGRAMS = unisolve
check_PROGRAMS = unisolve-test

unisolve_test_SOURCES = \
	unisolve-test.c \
	$(NULL)

unisolve_test_CFLAGS = \
	$(UNISOLVE_CFLAGS) \
	$(NULL)

unisolve_CFLAGS = \
	$(UNISOLVE_CFLAGS) \
	$(NULL)

unisolve_SOURCES = \
	unisolve.c \
	$(NULL)

unisolve_LDFLAGS = \
	$(UNISOLVE_LIBS)
	$(NULL)

unisolve_test_LDFLAGS = \
	$(UNISOLVE_LIBS) \
	$(NULL)

unisolve_LDADD = \
	$(top_builddir)/src/libmps/libmps.la \
	$(NULL)

unisolve_test_LDADD = \
	$(top_builddir)/src/libmps/libmps.la \
	$(NULL)
