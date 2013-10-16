if HAVE_UI

bin_PROGRAMS = $(QT_PROGRAMS)
EXTRA_PROGRAMS = xmpsolve$(EXEEXT)

man_MANS = ${builddir}/xmpsolve.1

# Install desktop file
desktopdir = $(datadir)/applications
desktop_DATA = xmpsolve.desktop

# Install icon for xmpsolve
icondir = $(datadir)/icons
icon_DATA = xmpsolve.png

xmpsolve_CPPFLAGS = \
	-I${top_srcdir}/include \
	-I${top_builddir} \
	-I${builddir} \
	-I${srcdir} \
	$(QT_CFLAGS)

if ! WINDOWS_BUILD
xmpsolve_CPPFLAGS += -fPIC
endif

xmpsolve_LDFLAGS = \
	$(QT_LDFLAGS)

xmpsolve_LDADD = \
	$(QT_LIBS) \
	${top_builddir}/src/libmps/libmps.la

%_moc.cpp: @srcdir@/%.h
	$(AM_V_GEN)
	$(AM_V_at) $(MOC) $< -o $@

ui_%.h: @srcdir@/%.ui
	$(AM_V_GEN) 
	$(AM_V_at) $(UIC) $< -o $@

%_qrc.cpp: @srcdir@/%.qrc
	$(AM_V_GEN)
	$(AM_V_at) $(RCC) $< -o $@

# Windows resources
if HAVE_WINDRES
.rc.o:
	$(AM_V_GEN)
	$(AM_V_at) $(WINDRES) $< -o $@

%.o : %.rc
	$(AM_V_GEN)
	$(AM_V_at) $(WINDRES) $< -o $@

resources = xmpsolve.rc

endif

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

ui_sfiles = \
	mainwindow.ui \
	polfileeditordialog.ui

qrc_files = \
	resources.qrc

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
qrc_gen_files = $(qrc_files:.qrc=_qrc.cpp)

xmpsolve_SOURCES = \
	$(other_source_files) \
	$(pure_source_files) \
	$(qrc_files) \
	$(resources)

nodist_xmpsolve_SOURCES = \
	$(moc_files) \
	$(qrc_gen_files)

CLEANFILES = \
	$(moc_files) \
	$(ui_files) \
	$(qrc_gen_files) 

EXTRA_DIST = $(icon_DATA) \
	$(headers) \
	$(desktop_DATA) \
	$(ui_sfiles) \
	xmpsolve.1

BUILT_SOURCES = $(moc_files) \
	$(ui_files) \
	$(qrc_gen_files) 

endif # HAVE_UI
