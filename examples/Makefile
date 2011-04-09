CFLAGS=-L../src/ -I../include/
LDFLAGS=-lgmp -lm -lmps

all: root_of_unity rand_poly

root_of_unity: root_of_unity.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o root_of_unity root_of_unity.c

rand_poly: root_of_unity.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o rand_poly rand_poly.c

clean:
	rm -f rand_poly root_of_unity
