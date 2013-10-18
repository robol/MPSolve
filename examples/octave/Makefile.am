if BUILD_OCT_FILE

mps_roots.o:
	$(MKOCTFILE) $(GMP_CFLAGS) -I${top_srcdir}/include -c -o mps_roots.o ${srcdir}/mps_roots.cc

noinst_PROGRAMS = $(OCTAVE_MODULES_PROGRAMS)
noinst_DATA = mpsolve.tar.gz

EXTRA_PROGRAMS = mps_roots.oct

EXTRA_DIST = \
	package \
	package/DESCRIPTION \
	package/COPYING \
	package/src

CLEANFILES = mpsolve.tar.gz

octave_share_prefix = ${DESTDIR}${datadir}/octave/packages
octave_lib_prefix =   ${DESTDIR}${libdir}/octave/packages

clean-local:
	rm -r ${builddir}/package

install-exec-hook:
	mkdir -p ${octave_share_prefix}
	mkdir -p ${octave_lib_prefix}
	octave -q -H --eval "pkg prefix ${octave_share_prefix} ${octave_lib_prefix} ; pkg install ${builddir}/mpsolve.tar.gz"

uninstall-local:
	octave -q -H --eval "pkg prefix ${octave_share_prefix} ${octave_lib_prefix} ; pkg uninstall mpsolve"

mpsolve.tar.gz: mps_roots.oct
	test -d package || (cp -r ${srcdir}/package/ ${builddir} && chmod u+w package -R)
	test -f package/COPYING || cp ${srcdir}/package/COPYING package/
	test -f package/DESCRIPTION || cp ${srcdir}/package/DESCRIPTION package/
	test -d package/src || mkdir package/src
	cp ${builddir}/mps_roots.oct ${builddir}/package/src/
	tar czf ${builddir}/mpsolve.tar.gz ${builddir}/package

mps_roots_oct_SOURCES = mps_roots.cc
mps_roots_oct_LINK = $(MKOCTFILE) -o $@ -L${top_builddir}/src/libmps/.libs -lmps 

endif
