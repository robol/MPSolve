if BUILD_OCT_FILE

mps_roots.o:
	$(MKOCTFILE) $(GMP_CFLAGS) -I${top_srcdir}/include -c -o mps_roots.o ${srcdir}/mps_roots.cc

noinst_PROGRAMS = mps_roots.oct
noinst_DATA = mpsolve.tar.gz

EXTRA_DIST = \
	package-source \
	package-source/DESCRIPTION \
	package-source/COPYING

CLEANFILES = mpsolve.tar.gz

octave_share_prefix = ${DESTDIR}${datadir}/octave/packages
octave_lib_prefix =   ${DESTDIR}${libdir}/octave/packages

clean-local:
	rm -rf ${builddir}/package

install-exec-hook:
	mkdir -p ${octave_share_prefix}
	mkdir -p ${octave_lib_prefix}
	octave -q -H --eval "pkg prefix ${octave_share_prefix} ${octave_lib_prefix}; \
	pkg global_list ${DESTDIR}/${datadir}/octave/octave_packages; \
	pkg install -verbose -nodeps \
	${builddir}/mpsolve.tar.gz"

uninstall-local:
	octave -q -H --eval "pkg prefix ${octave_share_prefix} ${octave_lib_prefix}; \
	pkg global_list ${DESTDIR}/${datadir}/octave/octave_packages; \
	pkg uninstall -verbose mpsolve"
	( rmdir ${DESTDIR}/${datadir}/octave/packages && rm -rf ${DESTDIR}/${datadir}/octave ) || \
	echo "Octave package dir is not empty, not removing it"

mpsolve.tar.gz: mps_roots.oct
	$(AM_V_GEN)
	$(AM_V_at) test -d package || (cp -r ${srcdir}/package-source ${builddir}/package && chmod u+w package -R)
	$(AM_V_at) test -f package/COPYING || cp ${srcdir}/package-source/COPYING package/
	$(AM_V_at) test -f package/DESCRIPTION || cp ${srcdir}/package-source/DESCRIPTION package/
	$(AM_V_at) test -d package/src || mkdir package/src
	$(AM_V_at) cp ${builddir}/mps_roots.oct ${builddir}/package/src/
	$(AM_V_at) tar czf ${builddir}/mpsolve.tar.gz ${builddir}/package

mps_roots_oct_SOURCES = mps_roots.cc
mps_roots_oct_LINK = $(MKOCTFILE) -o $@ -L${top_builddir}/src/libmps/.libs -lmps 

endif
