PREFIX=/usr
SRC=src
INSTALL=install

all: unisolve rursolve shared_libs

unisolve: shared_libs
	+$(MAKE) -C $(SRC) unisolve
	cp $(SRC)/unisolve .

rursolve: shared_libs unisolve
	+$(MAKE) -C $(SRC) rursolve
	cp $(SRC)/rursolve .

shared_libs:
	+$(MAKE) -C $(SRC) shared_libs

headers:
	mkdir -p $(PREFIX)/include/mps
	$(INSTALL) -m 644 include/mps/*.h $(PREFIX)/include/mps

install_libs: shared_libs
	install -m 644 $(SRC)/libmps.so $(PREFIX)/lib

install: unisolve rursolve install_libs headers
	install -m 755 unisolve $(PREFIX)/bin
	install -m 755 rursolve $(PREFIX)/bin

uninstall:
	rm -f $(PREFIX)/bin/unisolve
	rm -f $(PREFIX)/bin/rursolve
	rm -rf $(PREFIX)/include/mps
	rm -f $(PREFIX)/lib/libmps.so

check:
	make -C $(SRC) check
	make -C $(SRC) chckall	

clean:
	make -C $(SRC) clean
	rm -f unisolve rursolve

documentation: src/mps_*.c
	doxygen Doxyfile
