if BUILD_MATLAB_TOOLBOX

mps_roots.o:
	@$(MEX) -I${top_srcdir}/include -D__MPS_MATLAB_MODE \
	-c ${srcdir}/mps_roots.c > /dev/null 

noinst_PROGRAMS = $(MATLAB_MODULES_PROGRAMS)
noinst_DATA = mpsolve-toolbox.tar.gz

EXTRA_PROGRAMS = mps_roots

EXTRA_DIST = mps_roots.m

clean-local:
	rm -f mps_roots.mex*
	rm -rf mpsolve-toolbox/ 
	rm -f mpsolve-toolbox.tar.gz

mps_roots_SOURCES = mps_roots.c
mps_roots_LINK = $(MEX) -output $@ -L${top_builddir}/src/libmps/.libs \
	-lmps -lm -lgmp > /dev/null

mpsolve-toolbox.tar.gz: mps_roots mps_roots.m
	$(eval TEMP_DIR := $(shell mktemp -d tmp.XXXX))
	$(AM_V_GEN)
	$(AM_V_at) mkdir -p $(TEMP_DIR)/mpsolve-toolbox
	$(AM_V_at) cp mps_roots.mex* $(TEMP_DIR)/mpsolve-toolbox/ 
	$(AM_V_at) cp ${srcdir}/mps_roots.m $(TEMP_DIR)/mpsolve-toolbox/ 
	$(AM_V_at) cd $(TEMP_DIR) && tar czf ${builddir}/mpsolve-toolbox.tar.gz mpsolve-toolbox/
	$(AM_V_at) rm -rf $(TEMP_DIR)

endif
