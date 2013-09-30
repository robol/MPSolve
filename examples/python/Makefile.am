AM_CPPFLAGS = -I$(top_srcdir)/include $(PYTHON_INCLUDE)

EXTRA_PYTHON = mpsolve.la
EXTRA_DIST = mpsolve.pyx mpsolve.pxd
EXTRA_LTLIBRARIES = mpsolve.la

pyexec_LTLIBRARIES = $(PYTHON_MODULES)

mpsolve_la_SOURCES = $(builddir)/mpsolve.c
mpsolve_la_LDFLAGS = $(MPSOLVE_LDFLAGS) -module -avoid-version
mpsolve_la_LIBADD = $(top_builddir)/src/libmps/libmps.la

$(builddir)/mpsolve.c: $(srcdir)/mpsolve.pyx $(srcdir)/mpsolve.pxd
	$(CYTHON) -o $@ $(srcdir)/mpsolve.pyx -I$(top_srcdir)/include
