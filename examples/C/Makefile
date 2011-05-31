all: mandelbrot root_of_unity
	@echo "Compilation finished."
	@echo "You can now launch the examples if you have"
	@echo "installed MPSolve system-wide or settings  "
	@echo "the environment variable LD_LIBRARY_PATH   "
	@echo "to ../../src"
	@echo ""
	@echo "For example, try:"
	@echo "  LD_LIBRARY_PATH=../../src ./root_of_unity"


CFLAGS=-I../../include -DNOMPTEMP
LDFLAGS=-L../../src
CC=gcc

mandelbrot: mandelbrot.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o mandelbrot -lm -lmps -lgmp -lpthread mandelbrot.c

root_of_unity: root_of_unity.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o root_of_unity -lm -lmps -lgmp -lpthread root_of_unity.c

clean:
	rm root_of_unity mandelbrot
