bin_PROGRAMS = quadratic-solver

quadratic_solver_SOURCES = \
	quadratic-evaluation.c \
	quadratic-poly.c \
	quadratic-newton.c \
	quadratic-solver.c \
	quadratic-starting.c

quadratic_solver_CFLAGS = \
	-I${top_srcdir}/include \
	-I${top_builddir}/include

quadratic_solver_LDADD = \
	${top_builddir}/src/libmps/libmps.la \
	$(GMP_LIBS)

man_MANS = ${builddir}/quadratic-solver.1

EXTRA_DIST = \
	quadratic-poly.h 
