bin_PROGRAMS = $(QT4_PROGRAMS)
EXTRA_PROGRAMS = xmpsolve

xmpsolve_SOURCES = 

Makefile.qmake: xmpsolve.pro
	$(QMAKE_QT4) xmpsolve.pro

xmpsolve: Makefile.qmake
	+make -f Makefile.qmake xmpsolve	

clean-local:
	+make -f Makefile.qmake clean

distclean-local: clean-local
	rm -f Makefile.qmake
	rm -f xmpsolve.pro
	
