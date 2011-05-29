CFLAGS=-L../src/ -I../include/
MOD_CFLAGS=-I../include
LDFLAGS=-lgmp -lm -lmps -lpthread
MOD_LDFLAGS=-lgmp -lm -lpthread ../src/libmps.a ../src/libxt.a
MKOCTFILE=mkoctfile
MEX=mex
CC=gcc

all: root_of_unity mandelbrot

octave: mps_roots.oct
	@echo "The Octave module has been complex and has been named"
	@echo "mps_roots.oct. You can place it in the folder where "
	@echo "you usually launch octave or in the folder where octave"
	@echo "keeps all the .oct files, that usually is something"
	@echo "similar to /usr/lib/octave/<version>/oct/<architecture>/"

matlab: mps_roots.mex
	@echo "The MATLAB (TM) module has been compiled and has be named"
	@echo "mps_roots.mexglx. You can use it by placing in the path from"
	@echo "where MATLAB (TM) loads its MEX custom files, that is:"
	@echo " matlab-root/toolbox/local/"
	@echo "or you can simply launch matlab in the same folder that"
	@echo "mps_roots.mexglx"

root_of_unity: root_of_unity.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o root_of_unity root_of_unity.c

mandelbrot: mandelbrot.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o mandelbrot mandelbrot.c

mps_roots.oct: octave-module.cc
	$(MKOCTFILE) $(MOD_CFLAGS) -D__UNDEF_CPLUSPLUS -D__USE_BOOL_AS_BOOLEAN -o mps_roots.oct octave-module.cc $(MOD_LDFLAGS)

mps_roots.mex: matlab-module.c
	$(MEX) -output mps_roots -D__USE_BOOL_AS_BOOLEAN -L../src -I../include \
	$(MOD_LDFLAGS) matlab-module.c

secular: secular.c shared_libs
	$(CC) -o secular $(CFLAGS) $(LDFLAGS) secular.c

shared_libs: ../src/libmps.so
	$(MAKE) -C .. shared_libs

clean:
	rm -f mandelbrot root_of_unity mps_roots.oct *~ *.o *.mex
