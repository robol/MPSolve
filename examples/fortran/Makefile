F95=f95
CFLAGS=-I../../include
LDFLAGS=-L../../src -lmps -lm -lgmp -lpthread -static

roots_of_unity: mps.mod mps_impl.c roots_of_unity.f90
	$(F95) -o roots_of_unity roots_of_unity.f90 mps.f90 mps_impl.c \
	$(CFLAGS) $(LDFLAGS)

mps.mod: mps.f90
	$(F95) -c mps.f90

clean:
	rm -f *.mod roots_of_unity
