############################################################
##       Multiprecision Polynomial Solver (MPSolve)       ##
##                Version 2.2, May 2001                   ##
##                                                        ##
##                      Written by                        ##
##       Dario Andrea Bini and Giuseppe Fiorentino        ##
##       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        ##
##                                                        ##
## (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 ##
############################################################

# GNU compiler and linker invocation command
CC = gcc
LD = $(CC)

# path to GMP library and header file 
# (in case you want to use a GMP version different from installed)
GMPDIR = Gmp

# compiler flags
CFLAGS = -O2 -fPIC
# CFLAGS = -g -O0 -Wall -pedantic -fPIC
CPPFLAGS = -I$(GMPDIR)
#
# you can also use:
#  -DNOMPTEMP 	to disable the memory manager for MP temporary variables
#  -DRAND_VAL=x so that x will the return value for all rand() calls  
#

# linker flags
LDFLAGS = -L$(GMPDIR)
LOADLIBES = -lgmp -lm
#LOADLIBES = -lgmp -lm -lefence

# operating system dependent commands
AR = ar
ARFLAGS = rcs
RANLIB = ranlib
STRIP = strip
DIFF = diff
ECHO = echo
CAT = cat
RM = rm -f
SH = sh
DOXYGEN = doxygen

# operating system dependent standards
#HERE = 
HERE = ./
#EXE = .exe
EXE = 


# object files
LIBXTOBJ = tools.o mt.o gmptools.o mpc.o link.o mptemp.o
LIBMPSOBJ = mps_poly.o mps_data.o mps_sort.o mps_aber.o mps_clus.o mps_cnvx.o mps_impr.o mps_stio.o mps_newt.o mps_solv.o mps_star.o mps_test.o mps_touc.o mps_user.o mps_opts.o mps_main.o mps_defaults.o
RUROBJ = rur_horn.o rur_main.o


# main targets
unisolve: unisolve.o libmps.a libxt.a

rursolve: rursolve.o $(RUROBJ) libmps.a -lxt

shared_libs: libmps.so libxt.so

all: unisolve rursolve rurconv shared_libs

doc: 
	$(DOXYGEN) Doxyfile


# build libraries
libs: libmps.a libxt.a

libxt.a: $(LIBXTOBJ)
	$(AR) $(ARFLAGS) libxt.a $(LIBXTOBJ)
#	$(RANLIB) libxt.a

libmps.a: $(LIBMPSOBJ)
	$(AR) $(ARFLAGS) libmps.a $(LIBMPSOBJ)
#	$(RANLIB) libmps.a

libmps.so: $(LIBMPSOBJ)
	gcc -o libmps.so -shared -fPIC $(LIBMPSOBJ)

libxt.so: $(LIBXTOBJ)
	gcc -o libxt.so -shared -fPIC $(LIBXTOBJ)

install: libxt.so libmps.so
	install -m 644 libxt.so /usr/local/lib/
	install -m 644 libmps.so /usr/local/lib/
	install -m 644 mps.h /usr/local/include
	install -m 644 mps_poly.h /usr/local/include

# strip executables
strip:
	if [ -f unisolve$(EXE) ]; then $(STRIP) unisolve$(EXE); fi
	if [ -f rursolve$(EXE) ]; then $(STRIP) rursolve$(EXE); fi
	if [ -f rurconv$(EXE)  ]; then $(STRIP) rurconv$(EXE);  fi


# check
check: check1 check2 unisolve
	@$(ECHO) Good, unisolve seems to be ok!

check1: unisolve
	@$(ECHO) Test 1
	@$(ECHO) these roots:
	$(HERE)unisolve -Ga -Db -o100 Data/test.pol
	@$(ECHO) should be numerically equal to:
	$(CAT) Results/test.res
	@$(ECHO)

check2: unisolve
	@$(ECHO) Test 2
	$(HERE)unisolve -Ga -Db -o100 Data/wilk20.pol > Results/wilk20.chk
	$(DIFF) Results/wilk20.chk Results/wilk20.res
	@$(ECHO)

rurcheck: rursolve rurconv
	$(HERE)rurconv test.pl | $(HERE)rursolve


# cleanup
clean:
	$(RM) *.o stats timings roots testsuite *.bak *~ core *.res *.so

distclean: clean
	$(RM) libxt.a libmps.a unisolve$(EXE) rursolve$(EXE) rurconv$(EXE)

# test unisolve on all polynomials in the Data directory
# with a given target and set of options (see below)
testall: unisolve
	$(SH) Maketest $(TARGET) $(TESTOPT)

# test targets (see the maketest script for more details)
#TARGET = roots
#TARGET = time
TARGET = both
#TARGET = show
#TARGET = suite

# test options for unisolve (give your own or choose one below)
# TESTOPT1 has been used to compute the included results
TESTOPT = $(TESTOPT1)

# some sample options for benchmarks
TESTOPT1  = -Ga -o50
TESTOPT2  = -Gc -Si -o1000
TESTOPT3  = -Gc -SI -o1000
TESTOPT4  = -Gc -Si -o1000
TESTOPT5  = -Gc -So -o1000
TESTOPT5  = -Gc -Sl -o1000
TESTOPT6  = -Gc -Sr -o1000
TESTOPT7  = -Gc -Su -o1000
TESTOPT8  = -Gc -Sd -o1000
TESTOPT9  = -Gi -SR -o1000
TESTOPT10 = -Gi -SI -o1000
TESTOPT11 = -Gi -Si -o1000
TESTOPT12 = -Gi -So -o1000
TESTOPT13 = -Gi -Sl -o1000
TESTOPT14 = -Gi -Sr -o1000
TESTOPT15 = -Gi -Su -o1000
TESTOPT16 = -Gi -Sd -o1000
TESTOPT17 = -Gi -SR -Dr -o1000
TESTOPT18 = -Gi -SI -Di -o1000
TESTOPT19 = -Gi -Sa -Db -o1000
TESTOPT20 = -Gi -Sa -M+ -o1000

#
# ------------------------------
# targets for MPSolve developers
# ------------------------------
#
TOUCH = touch
DIFF = diff

# check targets
ress := $(wildcard Results/*.res)
chks := $(patsubst %.res, %.chk, $(ress))
pols := $(patsubst Results/%.res, Data/%.pol, $(ress))

setup: touchpol touchres

# run unisolve over all test polynomials with results
runall: $(pols)
	@$(ECHO) OK!

Data/%.pol: unisolve
	$(HERE)unisolve -d $(TESTOPT1) $@
	$(TOUCH) $@

touchpol:
	$(TOUCH) Data/*.pol

# (re)create all result files
results: $(ress)
	@$(ECHO) all results computed!

Results/%.res: Data/%.pol
#	$(HERE)unisolve $(TESTOPT1) $< > $@
	$(HERE)unisolve $(TESTOPT1) $< | tee $@

touchres:
	$(TOUCH) Results/*.res

# perform check against all test polynomials
chckall: $(chks)
	@$(ECHO) all tests succeeded!

Results/%.chk: Results/%.res
	$(HERE)unisolve $(TESTOPT1) Data/$*.pol > $@
	$(DIFF) $@ Results/$*.res

rmchk: touchpol touchres
	$(RM) Results/*.chk Results/*.gpt
	@$(ECHO) all check files removed

# version control and distribution targets 
CI = ci 
CO = co -q
COF = co -q -f
VERS = 2.2

co:
	$(CO) RCS/*,v ../XMT/RCS/*,v

ci: clean
	$(CI) RCS/*,v ../XMT/RCS/*,v

cof:
	$(COF) RCS/*,v ../XMT/RCS/*,v

version:
	ci -f -u$(VER) -m$(VERS) RCS/*,v ../XMT/RCS/*,v

getver:
	$(CO) -r$(VER) RCS/*,v ../XMT/RCS/*,v

# create distribution package
mantainerclean: distclean rmchk touchpol
	$(RM) Gmp/*

tgz: mantainerclean
	cd ..; \
	tar cvz --exclude RCS --exclude Tests --exclude Save \
	    --exclude Results/save -f MPSolve-$(VERS).tgz \
	    MPSolve-$(VERS)/

distrib: mantainerclean ci co tgz
