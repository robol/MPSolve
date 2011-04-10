CFLAGS=-L../src/ -I../include/
LDFLAGS=-lgmp -lm -lmps

all: root_of_unity mandelbrot

root_of_unity: root_of_unity.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o root_of_unity root_of_unity.c

mandelbrot: mandelbrot.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o mandelbrot mandelbrot.c

clean:
	rm -f mandelbrot root_of_unity
