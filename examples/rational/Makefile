RATIONAL_POLY_SOURCES = rational-poly.c \
	rational-poly-evaluation.c \
	rational-newton.c 

RATIONAL_POLY_OBJECTS = $(RATIONAL_POLY_SOURCES:.c=.o)

MEX=mkoctfile --mex

MATLAB_MEX=mex

# Please adjust the CFLAGS to your needs. This Makefile will build MEX files
# for Octave and MATLAB (if run with make matlab) but only when libmps and the
# related headers are available in the system. 
CFLAGS=-I../../include -I../../_build/include
rational_roots.mex: CFLAGS+=-D__MPS_NOT_DEFINE_BOOL

all: librational_poly.so rational_roots.mex

librational_poly.so: rational-poly.h $(RATIONAL_POLY_OBJECTS)
	$(CC) $(CFLAGS) -shared -o librational_poly.so $(RATIONAL_POLY_OBJECTS) -L. -lmps -Wl,-R.

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c -fPIC $<

rational_roots.mex: rational-roots.c librational_poly.so
	$(MEX) $(CFLAGS) -o $@ $< -lrational_poly -lmps -L. -Wl,-R.

matlab: rational-roots.c
	$(MATLAB_MEX) $(CFLAGS) -D__MPS_NOT_DEFINE_BOOL -I../../include -I../../_build/include -o rational_roots rational-roots.c -L. -lrational_poly -lmps -Wl,-R.

clean:
	rm -f *.o *.mex* librational_poly.so

check: rational_roots.mex
	octave -q test_rational.m
