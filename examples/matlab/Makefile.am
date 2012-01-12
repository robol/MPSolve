mps_roots.o:
	$(MEX) -I${top_srcdir}/include -D__MPS_MATLAB_MODE \
	-c ${srcdir}/mps_roots.c > /dev/null 

noinst_PROGRAMS = $(MATLAB_MODULES_PROGRAMS)

EXTRA_PROGRAMS = mps_roots
clean-local:
	rm -f mps_roots.mex*

mps_roots_SOURCES = mps_roots.c
mps_roots_LINK = $(MEX) -output $@ -L${top_builddir}/src/libmps/.libs \
	-lmps -lm -lgmp > /dev/null
