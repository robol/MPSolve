noinst_PROGRAMS = $(ROOTS_OF_UNITY_PROGRAMS)
EXTRA_PROGRAMS = roots_of_unity

roots_of_unity_SOURCES = roots_of_unity.c

roots_of_unity_CFLAGS = $(MPSOLVE_CFLAGS) \
	$(GTK_CFLAGS) \
	$(CAIRO_CFLAGS) \
	-I$(top_srcdir)/include

roots_of_unity_LDFLAGS = $(MPSOLVE_LDFLAGS) \
    $(CAIRO_LIBS) \
	$(GTK_LIBS) 

roots_of_unity_LDADD  = $(top_builddir)/src/libmps/libmps.la
