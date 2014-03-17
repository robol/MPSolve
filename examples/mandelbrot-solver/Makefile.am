bin_PROGRAMS = mandelbrot-solver

mandelbrot_solver_SOURCES = \
	mandelbrot-evaluation.c \
	mandelbrot-poly.c \
	mandelbrot-newton.c \
	mandelbrot-solver.c \
	mandelbrot-starting.c

mandelbrot_solver_CFLAGS = \
	-I${top_srcdir}/include \
	-I${top_builddir}/include

mandelbrot_solver_LDADD = \
	${top_builddir}/src/libmps/libmps.la \
	$(GMP_LIBS)

man_MANS = ${builddir}/mandelbrot-solver.1

EXTRA_DIST = \
	mandelbrot-poly.h 
