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

install_libs: shared_libs
	mkdir -p $(DESTDIR)/libs/mps
	$(INSTALL) -m 644 include/mps/*.h $(DESTDIR)/lib/mps

install: unisolve rursolve install_libs
	install -m 755 unisolve $(DESTDIR)/$(PREFIX)/bin
	install -m 755 rursolve $(DESTDIR)/$(PREFIX)/bin

clean:
	make -C $(SRC) clean
	rm -f unisolve rursolve

doc:
	make -C $(SRC) doc
