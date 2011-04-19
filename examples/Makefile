CFLAGS=-L../src/ -I../include/
LDFLAGS=-lgmp -lm -lmps
MKOCTFILE=mkoctfile

all: root_of_unity mandelbrot

root_of_unity: root_of_unity.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o root_of_unity root_of_unity.c

mandelbrot: mandelbrot.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o mandelbrot mandelbrot.c

mps_roots.oct: mps_roots.cc
	$(MKOCTFILE) mps_roots.cc $(LDFLAGS)

clean:
	rm -f mandelbrot root_of_unity
