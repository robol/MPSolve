bin_PROGRAMS = $(QT4_PROGRAMS)
EXTRA_PROGRAMS = xmpsolve$(EXEEXT)

man_MANS = ${srcdir}/xmpsolve.1

# Install desktop file
desktopdir = $(datadir)/applications
desktop_DATA = xmpsolve.desktop

# Install icon for xmpsolve
icondir = $(datadir)/icons
icon_DATA = xmpsolve.png

xmpsolve_CPPFLAGS = \
	-I${top_srcdir}/include \
	-I${top_builddir} \
	$(QT4_CFLAGS)

xmpsolve_LDADD = \
	$(QT4_LIBS) \
	${top_builddir}/src/libmps/libmps.la

EXTRA_DIST = $(icon_DATA) \
	$(desktop_DATA) \
	xmpsolve.1

%_moc.cpp: %.h
	$(MOC_QT4) $< -o $@

ui_%.h: %.ui
	$(UIC_QT4) $< -o $@

pure_source_files = \
	polynomialparser.cpp  \
	polynomialsolver.cpp  \
	mainwindow.cpp  \
	rootsmodel.cpp  \
	polfileeditordialog.cpp  \
	polsyntaxhighlighter.cpp  \
	mpsolveworker.cpp  \
	rootsrenderer.cpp  

other_source_files = \
	main.cpp \
	polynomial.cpp  \
	monomial.cpp  \
	root.cpp 

ui_files = \
	ui_mainwindow.h \
	ui_polfileeditordialog.h 

headers = \
	rootsrenderer.h  \
	rootsmodel.h  \
	polfileeditordialog.h  \
	monomial.h  \
	polynomial.h  \
	root.h  \
	polynomialparser.h  \
	polsyntaxhighlighter.h  \
	mainwindow.h  \
	polynomialsolver.h  \
	mpsolveworker.h

moc_files = $(pure_source_files:.cpp=_moc.cpp)

xmpsolve_SOURCES = \
	$(other_source_files) \
	$(pure_source_files) \
	$(moc_files)

clean-local:
	rm -f *_moc.cpp
	rm ui_*.h

BUILT_SOURCES = $(moc_files) \
	$(ui_files)
