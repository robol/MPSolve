TESTS = $(CHECK_PROGRAMS)
check_PROGRAMS = $(CHECK_PROGRAMS)
EXTRA_PROGRAMS = check_unisolve check_secsolve

check_secsolve_SOURCES = check_secsolve.c check_implementation.c
check_secsolve_CFLAGS = $(CHECK_CFLAGS) -I${top_builddir}/include
check_secsolve_LDFLAGS = -L${top_builddir}/src/libmps -lmps -lm -lgmp -lcheck

check_unisolve_SOURCES = check_unisolve.c check_implementation.c
check_unisolve_CFLAGS = $(CHECK_CFLAGS) -I${top_builddir}/include
check_unisolve_LDFLAGS = -L${top_builddir}/src/libmps -lmps -lm -lgmp -lcheck
