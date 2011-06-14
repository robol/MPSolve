NULL = 

lib_LTLIBRARIES = libmps.la
libmps_la_LDFLAGS = -version-info 3:0:0

AM_CPPFLAGS= \
	$(MPS_CFLAGS) \
	-pthread \
	$(NULL)

libmps_la_SOURCES = \
	 mps_user.c \
	 mps_opts.c \
	 mps_solv.c \
	 mps_main.c \
	 mps_data.c \
	 mps_clus.c \
	 mps_mptemp.c \
	 mps_stio.c \
	 mps_gmptools.c \
	 mps_defaults.c \
	 mps_sort.c \
	 mps_tools.c \
	 mps_poly.c \
	 mps_star.c \
	 mps_cnvx.c \
	 mps_mt.c \
	 mps_aber.c \
	 mps_impr.c \
	 mps_test.c \
	 mps_newt.c \
	 mps_mpc.c \
	 mps_link.c \
	 mps_touc.c \
	 mps_interface.c \
	 mps_secular.c \
	$(NULL)

libmps_la_LIBADD = \
	-lgmp \
	-lm \
	$(NULL)
	
