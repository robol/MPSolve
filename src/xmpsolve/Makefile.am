bin_PROGRAMS = $(QT4_PROGRAMS)
EXTRA_PROGRAMS = xmpsolve

xmpsolve_SOURCES = 

man_MANS = ${srcdir}/xmpsolve.1

EXTRA_DIST = xmpsolve.1

# Install desktop file
desktopdir = $(datadir)/applications
desktop_DATA = xmpsolve.desktop

# Install icon for xmpsolve
icondir = $(datadir)/icons
icon_DATA = xmpsolve.png

XMPSOLVE_QT_SOURCES = polynomialparser.cpp \
	rootsrenderer.h \
	polynomialsolver.cpp \
	mainwindow.cpp \
	polynomial.cpp \
	monomial.h \
	mainwindow.ui \
	polynomial.h \
	root.h \
	monomial.cpp \
	mpsolveworker.cpp \
	root.cpp \
	rootsrenderer.cpp \
	polynomialparser.h \
	mainwindow.h \
	polynomialsolver.h \
	mpsolveworker.h \
	main.cpp \
	mainwindow.ui

EXTRA_DIST = $(icon_DATA) \
	$(desktop_DATA) \
	$(XMPSOLVE_QT_SOURCES)

Makefile.qmake: xmpsolve.pro
	$(QMAKE_QT4) xmpsolve.pro

xmpsolve: Makefile.qmake
	+make -f Makefile.qmake xmpsolve	

clean-local:
	+make -f Makefile.qmake clean

distclean-local: clean-local
	rm -f Makefile.qmake
	rm -f ${srcdir}/xmpsolve.pro
