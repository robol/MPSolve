if BUILD_OCT_FILE

pkg_verbose = $(pkg_verbose_@AM_V@)
pkg_verbose_ = $(pkg_verbose_@AM_DEFAULT_V@)
pkg_verbose_0 = @echo "  MKOCT    $@";

noinst_PROGRAMS = mps_roots.oct mps_polyeig_impl.oct mps_det.oct
noinst_DATA = mpsolve.tar.gz

EXTRA_DIST = \
	DESCRIPTION \
	mps_polyeig.m \
	octave_support.h

CLEANFILES = mpsolve.tar.gz

octave_share_prefix = ${DESTDIR}${datadir}/octave/packages
octave_lib_prefix =   ${DESTDIR}${libdir}/octave/packages

clean-local:
	rm -rf ${builddir}/package

install-exec-hook:
	mkdir -p ${octave_share_prefix}
	mkdir -p ${octave_lib_prefix}
	LD_LIBRARY_PATH=${top_builddir}/src/libmps/.libs \ 
	    ${OCTAVE} -q -H --eval "pkg prefix ${octave_share_prefix} ${octave_lib_prefix}; \
	    pkg global_list ${DESTDIR}/${datadir}/octave/octave_packages; \
	    pkg install -verbose -nodeps \
	    ${builddir}/mpsolve.tar.gz"

uninstall-local:
	LD_LIBRARY_PATH=${top_builddir}/src/libmps/.libs \ 
	    ${OCTAVE} -q -H --eval "pkg prefix ${octave_share_prefix} ${octave_lib_prefix}; \
	    pkg global_list ${DESTDIR}/${datadir}/octave/octave_packages; \
	    pkg uninstall -verbose mpsolve"
	( rmdir ${DESTDIR}/${datadir}/octave/packages && rm -rf ${DESTDIR}/${datadir}/octave ) || \
	echo "Octave package dir is not empty, not removing it"

mpsolve.tar.gz: $(noinst_PROGRAMS) 
	$(AM_V_GEN)
	$(AM_V_at) test -d package || mkdir package
	$(AM_V_at) test -f package/COPYING || cp ${top_srcdir}/COPYING package/
	$(AM_V_at) test -f package/DESCRIPTION || cp ${srcdir}/DESCRIPTION package/
	$(AM_V_at) test -d package/src || mkdir package/src
	$(AM_V_at) cp ${builddir}/mps_roots.oct ${builddir}/package/src/
	$(AM_V_at) cp ${builddir}/mps_polyeig_impl.oct ${builddir}/package/src/
	$(AM_V_at) cp ${builddir}/mps_det.oct ${builddir}/package/src/
	$(AM_V_at) cp ${srcdir}/mps_polyeig.m ${builddir}/package/src/
	$(AM_V_at) tar czf ${builddir}/mpsolve.tar.gz ${builddir}/package

MKOCTFILE_LD = $(MKOCTFILE) -o $@ -L${top_builddir}/src/libmps/.libs -lmps

AM_CPPFLAGS = -fPIC 

mps_roots_oct_SOURCES = mps_roots.cc octave_support.cpp
mps_polyeig_impl_oct_SOURCES = mps_polyeig_impl.cc octave_support.cpp

mps_roots_oct_LINK = $(MKOCTFILE_LD)
mps_polyeig_impl_oct_LINK = $(MKOCTFILE_LD) 

mps_det_oct_SOURCES = mps_det.cc
mps_det_oct_LINK = $(MKOCTFILE_LD) 

mps_%.o: ${srcdir}/mps_%.cc 
	$(pkg_verbose)$(MKOCTFILE) $(GMP_CFLAGS) -I${top_srcdir}/include -I${top_builddir}/include -c -o $@ $<

endif
