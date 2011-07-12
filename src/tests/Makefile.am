TESTS = check_secsolve check_unisolve
check_PROGRAMS = check_secsolve check_unisolve

check_secsolve_SOURCES = check_secsolve.c check_implementation.c
check_secsolve_CFLAGS = @CHECK_CFLAGS@ -I${top_builddir}/include
check_secsolve_LDFLAGS = @CHECK_LIBS@ -L${top_builddir}/src/libmps -lmps -lm -lgmp

check_unisolve_SOURCES = check_unisolve.c check_implementation.c
check_unisolve_CFLAGS = @CHECK_CFLAGS@ -I${top_builddir}/include
check_unisolve_LDFLAGS = @CHECK_LIBS@ -L${top_builddir}/src/libmps -lmps -lm -lgmp
