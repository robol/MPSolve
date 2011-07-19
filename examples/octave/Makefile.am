mps_roots.o:
	$(MKOCTFILE) -I${top_srcdir}/include -DNOMPTEMP -DDISABLE_DEBUG \
	-D__UNDEF_CPLUSPLUS -D__USE_BOOL_AS_BOOLEAN \
	-c -o mps_roots.o ${srcdir}/mps_roots.cc

noinst_PROGRAMS = $(OCTAVE_MODULES_PROGRAMS)

EXTRA_PROGRAMS = mps_roots.oct
clean-local:
	rm -f mps_roots.oct

mps_roots_oct_SOURCES = mps_roots.cc
mps_roots_oct_LINK = $(MKOCTFILE) -o $@ -L${top_builddir}/src/libmps/.libs \
	-lmps -lm -lgmp
