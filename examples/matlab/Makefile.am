mps_roots.o:
	$(MEX) -I${top_builddir}/include -DNOMPTEMP -DDISABLE_DEBUG \
	-D__USE_BOOL_AS_BOOLEAN \
	-c mps_roots.c

noinst_PROGRAMS = $(MATLAB_MODULES_PROGRAMS)

EXTRA_PROGRAMS = mps_roots
clean-local:
	rm -f mps_roots.mex*

mps_roots_SOURCES = mps_roots.c
mps_roots_LINK = $(MEX) -output $@ -L${top_builddir}/src/libmps/.libs \
	-lmps -lm -lgmp
