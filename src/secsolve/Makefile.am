NULL = 

bin_PROGRAMS = secsolve
check_PROGRAMS = secsolve-test

secsolve_CFLAGS = \
	$(SECSOLVE_CFLAGS) \
	$(NULL)

secsolve_SOURCES = \
	secsolve.c \
	$(NULL)

secsolve_LDFLAGS = \
	$(SECSOLVE_LIBS)

secsolve_test_LDFLAGS = \
	$(SECSOLVE_LIBS)

secsolve_LDADD = \
	${top_builddir}/src/libmps/libmps.la

secsolve_test_CFLAGS = \
	$(SECSOLVE_CFLAGS)

secsolve_test_SOURCES = \
	secsolve-test.c

secsolve_test_LDADD = \
	${top_builddir}/src/libmps/libmps.la
