AM_CPPFLAGS = -I$(top_srcdir)/include $(PYTHON_INCLUDE)

EXTRA_PYTHON = _mpsolve.la
EXTRA_DIST = mpsolve.pyx mpsolve.pxd
EXTRA_LTLIBRARIES = _mpsolve.la

pyexec_LTLIBRARIES = $(PYTHON_MODULES)

_mpsolve_la_SOURCES = $(builddir)/mpsolve.c
_mpsolve_la_LDFLAGS = $(AM_LDFLAGS) $(MPSOLVE_LDFLAGS) -avoid-version -module
_mpsolve_la_LIBADD = $(top_builddir)/src/libmps/libmps.la

$(builddir)/mpsolve.c: $(srcdir)/mpsolve.pyx $(srcdir)/mpsolve.pxd
	@$(CYTHON) -o $@ $(srcdir)/mpsolve.pyx -I$(top_srcdir)/include
