MKOCTFILE=mkoctfile
CFLAGS=-I../../include
LDFLAGS=-L../../src -lgmp -lm -lpthread ../../src/libmps.a

all: mps_roots.oct
	@echo "The Octave module has been complex and has been named"
	@echo "mps_roots.oct. You can place it in the folder where "
	@echo "you usually launch octave or in the folder where octave"
	@echo "keeps all the .oct files, that usually is something"
	@echo "similar to /usr/lib/octave/<version>/oct/<architecture>/"

mps_roots.oct: octave-module.cc
	$(MKOCTFILE) $(CFLAGS) -D__UNDEF_CPLUSPLUS -D__USE_BOOL_AS_BOOLEAN -o mps_roots.oct octave-module.cc $(LDFLAGS)

clean:
	rm mps_roots.oct

