all: mandelbrot root_of_unity

CFLAGS=-I../../include -DNOMPTEMP -g -O0
LDFLAGS=-lmps -lm -lgmp -lpthread -L../../src/libmps/.libs -static
CC=gcc

mandelbrot: mandelbrot.c
	$(CC) -o mandelbrot mandelbrot.c $(CFLAGS) $(LDFLAGS)


root_of_unity: root_of_unity.c
	$(CC) -o root_of_unity root_of_unity.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f root_of_unity mandelbrot
