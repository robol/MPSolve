if CHECK

check_PROGRAMS = check_convex check_context check_mpc check_matrix check_dpe \
	check_formal \
	check_multithread check_cluster check_chebyshev check_parser check_utils \
	check_monomial_poly check_list check_secsolve check_unisolve check_version

TESTS = $(check_PROGRAMS)  

AM_TEST_ENVIRONMENT = \
  SRCDIR=$(abs_srcdir) BUILDDIR=$(abs_builddir)

COMMON_LIBS = -lm $(GMP_LIBS) $(CHECK_LIBS) $(PTHREAD_LIBS) -lgmpxx
COMMON_CFLAGS = $(CHECK_CFLAGS) $(GMP_CFLAGS) $(PTHREAD_CFLAGS) -I${top_builddir}/include -I${top_srcdir}/include -D_MPS_PRIVATE -g -O0
COMMON_SOURCES = check_implementation.c check_implementation.h dummy.cpp
COMMON_LDADD = ${top_builddir}/src/libmps/libmpsprivate.la

 check_formal_SOURCES = check_formal.cpp $(COMMON_SOURCES)
 check_formal_CFLAGS = $(COMMON_CFLAGS)
 check_formal_CPPFLAGS = $(COMMON_CFLAGS)
 check_formal_LDFLAGS = $(COMMON_LIBS) 
 check_formal_LDADD = $(COMMON_LDADD)

 check_list_SOURCES = check_list.c $(COMMON_SOURCES)
 check_list_CFLAGS = $(COMMON_CFLAGS)
 check_list_LDFLAGS = $(COMMON_LIBS) 
 check_list_LDADD = $(COMMON_LDADD)

 check_convex_SOURCES = check_convex.c $(COMMON_SOURCES)
 check_convex_CFLAGS = $(COMMON_CFLAGS)
 check_convex_LDFLAGS = $(COMMON_LIBS) 
 check_convex_LDADD = $(COMMON_LDADD)

 check_parser_SOURCES = check_parser.c $(COMMON_SOURCES)
 check_parser_CFLAGS = $(COMMON_CFLAGS)
 check_parser_LDFLAGS = $(COMMON_LIBS) 
 check_parser_LDADD = $(COMMON_LDADD)

 check_monomial_poly_SOURCES = check_monomial_poly.c $(COMMON_SOURCES)
 check_monomial_poly_CFLAGS = $(COMMON_CFLAGS)
 check_monomial_poly_LDFLAGS = $(COMMON_LIBS) 
 check_monomial_poly_LDADD = $(COMMON_LDADD)

 check_utils_SOURCES = check_utils.c $(COMMON_SOURCES)
 check_utils_CFLAGS = $(COMMON_CFLAGS)
 check_utils_LDFLAGS = $(COMMON_LIBS) 
 check_utils_LDADD = $(COMMON_LDADD)

 check_mpc_SOURCES = check_mpc.c $(COMMON_SOURCES)
 check_mpc_CFLAGS = $(COMMON_CFLAGS)
 check_mpc_LDFLAGS = $(COMMON_LIBS) 
 check_mpc_LDADD = $(COMMON_LDADD)

 check_context_SOURCES = check_context.c $(COMMON_SOURCES)
 check_context_CFLAGS = $(COMMON_CFLAGS)
 check_context_LDFLAGS = $(COMMON_LIBS) 
 check_context_LDADD = $(COMMON_LDADD)

 check_chebyshev_SOURCES = check_chebyshev.c $(COMMON_SOURCES)
 check_chebyshev_CFLAGS = $(COMMON_CFLAGS)
 check_chebyshev_LDFLAGS = $(COMMON_LIBS)
 check_chebyshev_LDADD = $(COMMON_LDADD)

 check_matrix_SOURCES = check_matrix.c $(COMMON_SOURCES) 
 check_matrix_CFLAGS = $(COMMON_CFLAGS)
 check_matrix_LDFLAGS = $(COMMON_LIBS) 
 check_matrix_LDADD = $(COMMON_LDADD)

 check_secsolve_SOURCES = check_secsolve.c $(COMMON_SOURCES)
 check_secsolve_CFLAGS = $(COMMON_CFLAGS)
 check_secsolve_LDFLAGS = $(COMMON_LIBS)
 check_secsolve_LDADD = $(COMMON_LDADD)

 check_unisolve_SOURCES = check_unisolve.c $(COMMON_SOURCES)
 check_unisolve_CFLAGS = $(COMMON_CFLAGS)
 check_unisolve_LDFLAGS = $(COMMON_LIBS)
 check_unisolve_LDADD = $(COMMON_LDADD)

 check_dpe_SOURCES = check_dpe.c $(COMMON_SOURCES)
 check_dpe_CFLAGS = $(COMMON_CFLAGS)
 check_dpe_LDFLAGS = $(COMMON_LIBS)
 check_dpe_LDADD = $(COMMON_LDADD)

 check_multithread_SOURCES = check_multithread.c $(COMMON_SOURCES)
 check_multithread_CFLAGS = $(COMMON_CFLAGS)
 check_multithread_LDFLAGS = $(COMMON_LIBS)
 check_multithread_LDADD = $(COMMON_LDADD)

 check_cluster_SOURCES = check_cluster.c $(COMMON_SOURCES)
 check_cluster_CFLAGS = $(COMMON_CFLAGS)
 check_cluster_LDFLAGS = $(COMMON_LIBS)
 check_cluster_LDADD = $(COMMON_LDADD)

 check_version_SOURCES = check_version.c $(COMMON_SOURCES)
 check_version_CFLAGS = $(COMMON_CFLAGS)
 check_version_LDFLAGS = $(COMMON_LIBS)
 check_version_LDADD = $(COMMON_LDADD)

endif

EXTRA_DIST = \
	secsolve/rand15.pol \
	secsolve/test1000.pol \
	secsolve/test100.pol \
	secsolve/rand120.pol \
	secsolve/integer.pol \
	secsolve/simple.pol \
	secsolve/wilk20.pol \
	secsolve/wilk40.pol \
	secsolve/wilk80.pol \
	unisolve/lsr_24.pol \
	unisolve/toep1_128.pol \
	unisolve/kir1_20_mod.pol \
	unisolve/easy400.pol \
	unisolve/laguerre20.pol \
	unisolve/nroots800.pol \
	unisolve/lar3.pol \
	unisolve/geom2_20.pol \
	unisolve/sendra20.pol \
	unisolve/easy3200.pol \
	unisolve/laguerre40.pol \
	unisolve/nroots3200.pol \
	unisolve/mand255.pol \
	unisolve/hermite20.pol \
	unisolve/kir1_20.pol \
	unisolve/lar1.pol \
	unisolve/spiral25.pol \
	unisolve/chebyshev20.pol \
	unisolve/nrooti200.pol \
	unisolve/easy100.pol \
	unisolve/mand31.pol \
	unisolve/kam4.pol \
	unisolve/wilk20.pol \
	unisolve/wilk40.pol \
	unisolve/wilk80.pol \
	unisolve/kam1_1.pol \
	unisolve/chebyshev80.pol \
	unisolve/kam3_3.pol \
	unisolve/mult1.pol \
	unisolve/geom4_80.pol \
	unisolve/chebyshev40.pol \
	unisolve/legendre80.pol \
	unisolve/sendra320.pol \
	unisolve/lsr3.pol \
	unisolve/mig1_100_1.pol \
	unisolve/lsr4_3.pol \
	unisolve/lar2.pol \
	unisolve/hermite320.pol \
	unisolve/kir1_symb.pol \
	unisolve/mig1_50_1.pol \
	unisolve/geom4_40.pol \
	unisolve/nroots400.pol \
	unisolve/chrmc_d11.pol \
	unisolve/toep2_128.pol \
	unisolve/hermite160.pol \
	unisolve/chebyshev320.pol \
	unisolve/geom3_20.pol \
	unisolve/geom1_15.pol \
	unisolve/geom2_15.pol \
	unisolve/nroots200.pol \
	unisolve/umand255.pol \
	unisolve/kam2_1.pol \
	unisolve/easy800.pol \
	unisolve/chrma_d340.pol \
	unisolve/chrma22.pol \
	unisolve/kir1_40_mod.pol \
	unisolve/mand511.pol \
	unisolve/nrooti100.pol \
	unisolve/sparse6400.pol \
	unisolve/easy200.pol \
	unisolve/lsr1.pol \
	unisolve/kam2_3.pol \
	unisolve/mig1_20.pol \
	unisolve/curz40.pol \
	unisolve/exp200.pol \
	unisolve/nrooti1600.pol \
	unisolve/test.pol \
	unisolve/sparse100.pol \
	unisolve/mand1023.pol \
	unisolve/exp400.pol \
	unisolve/legendre320.pol \
	unisolve/laguerre320.pol \
	unisolve/nektarios.pol \
	unisolve/chebyshev160.pol \
	unisolve/mig1_200.pol \
	unisolve/kam1_3.pol \
	unisolve/kir1_40.pol \
	unisolve/chrmc_d43.pol \
	unisolve/geom1_10.pol \
	unisolve/mult4.pol \
	unisolve/lsr2.pol \
	unisolve/kam2_2.pol \
	unisolve/chrmc_d683.pol \
	unisolve/mand63.pol \
	unisolve/hermite40.pol \
	unisolve/umand2047.pol \
	unisolve/geom1_20.pol \
	unisolve/chrmc343.pol \
	unisolve/chrmc23.pol \
	unisolve/chrma_d20.pol \
	unisolve/kam1_2.pol \
	unisolve/umand127.pol \
	unisolve/easy1600.pol \
	unisolve/sendra160.pol \
	unisolve/wilk320.pol \
	unisolve/curz20.pol \
	unisolve/toep1_256.pol \
	unisolve/geom1_40.pol \
	unisolve/wilk160.pol \
	unisolve/mult2.pol \
	unisolve/nrooti3200.pol \
	unisolve/lsr_200.pol \
	unisolve/lsr4_2.pol \
	unisolve/wilk_mod.pol \
	unisolve/sparse200.pol \
	unisolve/sendra40.pol \
	unisolve/legendre40.pol \
	unisolve/curz160.pol \
	unisolve/chrmc_d171.pol \
	unisolve/geom3_10.pol \
	unisolve/chrma342.pol \
	unisolve/chrma86.pol \
	unisolve/spiral10.pol \
	unisolve/curz80.pol \
	unisolve/sparse800.pol \
	unisolve/nrooti6400.pol \
	unisolve/nrooti50.pol \
	unisolve/kir1_10.pol \
	unisolve/spiral15.pol \
	unisolve/mand127.pol \
	unisolve/hermite80.pol \
	unisolve/nroots50.pol \
	unisolve/umand1023.pol \
	unisolve/nrooti400.pol \
	unisolve/geom3_40.pol \
	unisolve/legendre160.pol \
	unisolve/lar4.pol \
	unisolve/geom3_80.pol \
	unisolve/nroots1600.pol \
	unisolve/mig1_200_1.pol \
	unisolve/spiral30.pol \
	unisolve/wilk80.pol \
	unisolve/trv_m.pol \
	unisolve/geom2_40.pol \
	unisolve/laguerre80.pol \
	unisolve/sparse400.pol \
	unisolve/geom2_10.pol \
	unisolve/umand31.pol \
	unisolve/sendra80.pol \
	unisolve/kam3_2.pol \
	unisolve/lsr4_1.pol \
	unisolve/sparse3200.pol \
	unisolve/geom4_10.pol \
	unisolve/mig1_500_1.pol \
	unisolve/legendre20.pol \
	unisolve/toep2_256.pol \
	unisolve/chrma_d84.pol \
	unisolve/nroots6400.pol \
	unisolve/lar1_200.pol \
	unisolve/mult3.pol \
	unisolve/exp50.pol \
	unisolve/kats8.pol \
	unisolve/lar5.pol \
	unisolve/nroots100.pol \
	unisolve/laguerre160.pol \
	unisolve/exp100.pol \
	unisolve/kam3_1.pol \
	unisolve/mig1_500.pol \
	unisolve/geom4_20.pol \
	unisolve/umand511.pol \
	unisolve/kir1_10_mod.pol \
	unisolve/sparse1600.pol \
	unisolve/wilk40.pol \
	unisolve/spiral20.pol \
	unisolve/nrooti800.pol \
	unisolve/umand63.pol \
	unisolve/mig1_100.pol \
	unisolve/demi20.pol \
	../results/secsolve/test1000.res \
	../results/secsolve/rand15.res \
	../results/secsolve/test100.res \
	../results/secsolve/rand120.res \
	../results/secsolve/integer.res \
	../results/secsolve/simple.res \
	../results/secsolve/wilk20.res \
	../results/secsolve/wilk40.res \
	../results/secsolve/wilk80.res \
	../results/unisolve/mand63.res \
	../results/unisolve/kam2_1.res \
	../results/unisolve/kir1_10.res \
	../results/unisolve/kir1_20.res \
	../results/unisolve/kir1_40.res \
	../results/unisolve/spiral10.res \
	../results/unisolve/spiral20.res \
	../results/unisolve/nroots50.res \
	../results/unisolve/wilk40.res \
	../results/unisolve/wilk80.res \
	../results/unisolve/lar3.res \
	../results/unisolve/trv_m.res \
	../results/unisolve/lar2.res \
	../results/unisolve/toep1_128.res \
	../results/unisolve/mult1.res \
	../results/unisolve/exp50.res \
	../results/unisolve/mand127.res \
	../results/unisolve/kam3_2.res \
	../results/unisolve/kam3_3.res \
	../results/unisolve/umand31.res \
	../results/unisolve/kam1_2.res \
	../results/unisolve/kam4.res \
	../results/unisolve/lar1.res \
	../results/unisolve/spiral20.res \
	../results/unisolve/wilk20.res \
	../results/unisolve/test.res \
	../results/unisolve/lar1_200.res \
	../results/unisolve/kam1_3.res \
	../results/unisolve/kam3_1.res \
	../results/unisolve/kam2_2.res \
	../results/unisolve/exp100.res \
	../results/unisolve/kam1_1.res \
	../results/unisolve/kam2_3.res \
	../results/unisolve/mig1_100.res \
	../results/unisolve/mig1_200.res \
	../results/unisolve/lsr_24.res \
	../results/unisolve/mig1_500_1.res \
	../results/unisolve/demi20.res
