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
LDFLAGS=-lmps -lm -lgmp -lpthread -L../../src -static
CC=gcc

mandelbrot: mandelbrot.c
	$(CC) -o mandelbrot mandelbrot.c $(CFLAGS) $(LDFLAGS)


root_of_unity: root_of_unity.c
	$(CC) -o root_of_unity root_of_unity.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f root_of_unity mandelbrot
