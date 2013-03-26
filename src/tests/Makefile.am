TESTS = $(CHECK_PROGRAMS)
check_PROGRAMS = $(CHECK_PROGRAMS)
EXTRA_PROGRAMS = check_secsolve check_unisolve check_dpe test_mpsolve check_multithread check_cluster

 check_secsolve_SOURCES = check_secsolve.c check_implementation.c check_implementation.h
 check_secsolve_CFLAGS = $(CHECK_CFLAGS) $(GMP_CFLAGS) $(PTHREAD_CFLAGS) -I${top_srcdir}/include
 check_secsolve_LDFLAGS = ${top_builddir}/src/libmps/libmps.la -lm $(GMP_LIBS) $(CHECK_LIBS) $(PTHREAD_LIBS)

 check_unisolve_SOURCES = check_unisolve.c check_implementation.c check_implementation.h
 check_unisolve_CFLAGS = $(CHECK_CFLAGS) $(GMP_CFLAGS) $(PTHREAD_CFLAGS) -I${top_srcdir}/include
 check_unisolve_LDFLAGS = ${top_builddir}/src/libmps/libmps.la -lm $(GMP_LIBS) $(CHECK_LIBS) $(PTHREAD_LIBS)

 check_dpe_SOURCES = check_dpe.c check_implementation.c check_implementation.h
 check_dpe_CFLAGS = $(CHECK_CFLAGS) $(GMP_CFLAGS) $(PTHREAD_CFLAGS) -I${top_srcdir}/include
 check_dpe_LDFLAGS = ${top_builddir}/src/libmps/libmps.la -lm $(GMP_LIBS) $(CHECK_LIBS) $(PTHREAD_LIBS)

 test_mpsolve_SOURCES = test_mpsolve.c check_implementation.h check_implementation.c
 test_mpsolve_CFLAGS = ${CHECK_CFLAGS} $(GMP_CFLAGS) $(PTHREAD_CFLAGS) -I${top_srcdir}/include
 test_mpsolve_LDFLAGS = ${top_builddir}/src/libmps/libmps.la -lm $(GMP_LIBS) $(PTHREAD_LIBS)

 check_multithread_SOURCES = check_multithread.c
 check_multithread_CFLAGS = $(CHECK_CFLAGS) $(GMP_CFLAGS) $(PTHREAD_CFLAGS) -I${top_srcdir}/include -D_MPS_PRIVATE
 check_multithread_LDFLAGS = ${top_builddir}/src/libmps/libmps.la -lm $(GMP_LIBS) $(CHECK_LIBS) $(PTHREAD_LIBS)

 check_cluster_SOURCES = check_cluster.c check_implementation.c
 check_cluster_CFLAGS = $(CHECK_CFLAGS) $(GMP_CFLAGS) $(PTHREAD_CFLAGS) -I${top_srcdir}/include -D_MPS_PRIVATE
 check_cluster_LDFLAGS = ${top_builddir}/src/libmps/libmps.la -lm $(GMP_LIBS) $(CHECK_LIBS) $(PTHREAD_LIBS)

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
	../results/unisolve/mig1_500_1.res 
