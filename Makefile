PREFIX=/usr
DESTDIR=
SRC=src
INSTALL=install

all: unisolve rursolve shared_libs

unisolve:
	make -C $(SRC) unisolve
	cp $(SRC)/unisolve .

rursolve:
	make -C $(SRC) rursolve
	cp $(SRC)/rursolve .

shared_libs:
	make -C $(SRC) shared_libs

headers:
	mkdir -p $(DESTDIR)/$(PREFIX)/include/mps
	$(INSTALL) -m 644 include/mps/*.h $(DESTDIR)/$(PREFIX)/include/mps

install_libs: shared_libs
	install -m 644 $(SRC)/libmps.so $(DESTDIR)/$(PREFIX)/lib
	install -m 644 $(SRC)/libxt.so $(DESTDIR)/$(PREFIX)/lib

install: unisolve rursolve install_libs headers
	install -m 755 unisolve $(DESTDIR)/$(PREFIX)/bin
	install -m 755 rursolve $(DESTDIR)/$(PREFIX)/bin

uninstall:
	rm -f $(DESTDIR)/$(PREFIX)/bin/unisolve
	rm -f $(DESTDIR)/$(PREFIX)/bin/rursolve
	rm -rf $(DESTDIR)/$(PREFIX)/include/mps
	rm -f $(DESTDIR)/$(PREFIX)/lib/libmps.so
	rm -f $(DESTDIR)/$(PREFIX)/lib/libxt.so

check:
	make -C $(SRC) check
	make -C $(SRC) chckall	

clean:
	make -C $(SRC) clean
	rm -f unisolve rursolve

doc:
	make -C $(SRC) doc
