CC=gcc

roots_of_unity: roots_of_unity.c
	$(CC) -o roots_of_unity roots_of_unity.c -L../../src/libmps/.libs -lmps `pkg-config --cflags --libs gtk+-3.0` -I../../include

clean:
	rm -f roots_of_unity
