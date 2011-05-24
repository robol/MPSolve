PREFIX=/usr
SRC=src
INSTALL=install
CC=gcc
LD=gcc
DESTDIR=
LDFLAGS=-lgmp -lm -lpthread
MAKEFLAGS=-j

# CFLAGS for the file in src/. 
# The debug ones are used if the debug variable
# is set, otherwise normal complation is performed. 
# You can also use:
#  -DNOMPTEMP   to disable the memory manager for MP temporary variables
#  -DRAND_VAL=x so that x will the return value for all rand() calls

# Default CFLAGS
CFLAGS=-O2 -ffast-math -fPIC -I../include -std=c99 -DNOMPTEMP

# Set CFLAGS for specific targets
debug: CFLAGS=-O0 -g -fPIC -I../include -Wall -std=c99 -DNOMPTEMP
release: CFLAGS=-O2 -ffast-math -fPIC -DDISABLE_DEBUG -I../include -std=c99 -DNOMPTEMP

# Export CFLAGS and CC for the submakes. 
export CFLAGS
export LDFLAGS
export CC

#
# Targets
# 
all: unisolve rursolve secsolve shared_libs

debug: all

release: all

unisolve: shared_libs
	+$(MAKE) -C $(SRC) unisolve
	cp $(SRC)/unisolve .

rursolve: shared_libs unisolve
	+$(MAKE) -C $(SRC) rursolve
	cp $(SRC)/rursolve .

secsolve: shared_libs rursolve
	+$(MAKE) -C $(SRC) secsolve
	cp $(SRC)/secsolve .

shared_libs:
	+$(MAKE) -C $(SRC) shared_libs

headers:
	mkdir -p $(PREFIX)/include/mps
	$(INSTALL) -m 644 include/mps/*.h $(DESTDIR)$(PREFIX)/include/mps

install_libs: shared_libs
	install -m 644 $(SRC)/libmps.so $(DESTDIR)$(PREFIX)/lib

install: unisolve rursolve install_libs headers
	install -m 755 unisolve $(DESTDIR)$(PREFIX)/bin
	install -m 755 rursolve $(DESTDIR)$(PREFIX)/bin

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/unisolve
	rm -f $(DESTDIR)$(PREFIX)/bin/rursolve
	rm -rf $(DESTDIR)$(PREFIX)/include/mps
	rm -f $(DESTDIR)$(PREFIX)/lib/libmps.so

check:
	make -C $(SRC) check
	make -C $(SRC) chckall	

clean:
	make -C $(SRC) clean
	rm -f unisolve rursolve secsolve

distclean: clean
	rm -rf doc/html doc/latex

documentation: src/mps_*.c
	doxygen Doxyfile
