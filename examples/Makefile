CFLAGS=-L../src/ -I../include/
LDFLAGS=-lgmp -lm -lmps
MKOCTFILE=mkoctfile

all: root_of_unity mandelbrot

root_of_unity: root_of_unity.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o root_of_unity root_of_unity.c

mandelbrot: mandelbrot.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o mandelbrot mandelbrot.c

mps_roots.oct: octave-module.cc
	$(MKOCTFILE) -o mps_roots.oct octave-module.cc $(LDFLAGS)

clean:
	rm -f mandelbrot root_of_unity mps_roots.oct *~ *.o
