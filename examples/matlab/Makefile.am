if BUILD_MATLAB_TOOLBOX

%.o: %.c
	@$(MEX) $(MEXOPTS) $(GMP_CFLAGS) -I${top_srcdir}/include -g \
	-I${top_builddir}/include \
	-c $< > /dev/null

noinst_PROGRAMS = $(MATLAB_MODULES_PROGRAMS)
noinst_DATA = mpsolve-toolbox.tar.gz

EXTRA_PROGRAMS = mps_roots_double mps_roots_string mps_kostlan mps_interp

EXTRA_DIST = mps_roots.m load_mpsolve.m

clean-local:
	rm -f mps_roots_double.mex* mps_roots_string.mex* mps_koslan.mex* mps_interp.mex*
	rm -rf mpsolve-toolbox/ 
	rm -f mpsolve-toolbox.tar.gz

MPS_LINK_COMMAND = $(MEX) $(MEXOPTS) -output $@ -L${top_builddir}/src/libmps/.libs \
	$(GMP_LIBS) -lmps -lm -lgmp -Wl,-R. > /dev/null

mps_kostlan_SOURCES = mps_kostlan.c mps_option_parser.c mps_option_parser.h
mps_kostlan_LINK = ${MPS_LINK_COMMAND}

mps_interp_SOURCES = mps_interp.c mps_option_parser.c mps_option_parser.h
mps_interp_LINK = ${MPS_LINK_COMMAND} -lmpfr

mps_roots_double_SOURCES = mps_roots_double.c mps_option_parser.c mps_option_parser.h
mps_roots_double_LINK = ${MPS_LINK_COMMAND} 

mps_roots_string_SOURCES = mps_roots_string.c mps_option_parser.c mps_option_parser.h
mps_roots_string_LINK = ${MPS_LINK_COMMAND}

mpsolve-toolbox.tar.gz: mps_roots_double mps_roots_string mps_kostlan mps_interp mps_roots.m load_mpsolve.m ${top_builddir}/src/libmps/.libs/libmps.so.3
	$(eval TEMP_DIR := $(shell mktemp -d tmp.XXXX))
	$(AM_V_GEN)
	$(AM_V_at) mkdir -p $(TEMP_DIR)/mpsolve-toolbox/private
	$(AM_V_at) cp mps_roots*.mex* $(TEMP_DIR)/mpsolve-toolbox/private/
	$(AM_V_at) cp ${top_builddir}/src/libmps/.libs/libmps.so.3 $(TEMP_DIR)/
	$(AM_V_at) cp ${srcdir}/mps_roots.m mps_kostlan.mex* mps_interp.mex* $(TEMP_DIR)/mpsolve-toolbox/ 
	$(AM_V_at) cp ${srcdir}/load_mpsolve.m $(TEMP_DIR)/
	$(AM_V_at) cd $(TEMP_DIR) && tar czf ${builddir}/mpsolve-toolbox.tar.gz *
	$(AM_V_at) mv $(TEMP_DIR)/mpsolve-toolbox.tar.gz ${builddir}/
	$(AM_V_at) rm -rf $(TEMP_DIR)

endif
