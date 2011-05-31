F95=f95
CFLAGS=-I../../include
LDFLAGS=-L../../src

roots_of_unity:
	$(F95) -o roots_of_unity roots_of_unity.f90 mps.f90 mps_impl.c -lm -lmps -lgmp -lpthread \
	$(CFLAGS) $(LDFLAGS)

clean:
	rm *.mod roots_of_unity
