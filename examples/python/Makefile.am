AM_CPPFLAGS = -I$(top_srcdir)/include $(PYTHON_INCLUDE)

EXTRA_DIST = mpsolve.pyx mpsolve.pxd

if HAVE_CYTHON
  pyexec_LTLIBRARIES = mpsolve.la
else
  pyexec_LTLIBRARIES =
endif

nodist_mpsolve_la_SOURCES = $(builddir)/mpsolve.c
mpsolve_la_LDFLAGS = $(AM_LDFLAGS) $(MPSOLVE_LDFLAGS) -avoid-version -module
mpsolve_la_LIBADD = $(top_builddir)/src/libmps/libmps.la

if HAVE_CYTHON
$(builddir)/mpsolve.c: $(srcdir)/mpsolve.pyx $(srcdir)/mpsolve.pxd
	@$(CYTHON) -o $@ $(srcdir)/mpsolve.pyx -I$(top_srcdir)/include
endif
