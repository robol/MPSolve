MEX=mex -g
CFLAGS=-D__USE_BOOL_AS_BOOLEAN -DNOMPTEMP -I../../include
LDFLAGS=-L../../src/libmps/.libs -lm -lgmp -lpthread -lmps

all: matlab
	@echo "The MATLAB (TM) module has been compiled and has be named"
	@echo "mps_roots.mexglx. You can use it by placing in the path from"
	@echo "where MATLAB (TM) loads its MEX custom files, that is:"
	@echo " matlab-root/toolbox/local/"
	@echo "or you can simply launch matlab in the same folder that"
	@echo "mps_roots.mexglx"

matlab:
	$(MEX) -output mps_roots matlab-module.c $(LDFLAGS) $(CFLAGS)

clean:
	rm -f mps_roots*
