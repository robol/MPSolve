if HAVE_UI

bin_PROGRAMS = $(QT_PROGRAMS)
EXTRA_PROGRAMS = xmpsolve

man_MANS = ${builddir}/xmpsolve.1

# Install desktop file
desktopdir = $(datadir)/applications
desktop_DATA = xmpsolve.desktop

# Install icon for xmpsolve
icondir = $(datadir)/icons
icon_DATA = xmpsolve.png

xmpsolve_CPPFLAGS = \
	-I${top_srcdir}/include \
	-I${top_builddir}/include \
	-I${builddir} \
        -I${top_builddir}/include \
	-I${srcdir} \
	$(QT_CFLAGS)

# if ! WINDOWS_BUILD
xmpsolve_CPPFLAGS += -fPIC
# endif

xmpsolve_LDFLAGS = \
	$(QT_LDFLAGS)

xmpsolve_LDADD = \
	$(QT_LIBS) \
	${top_builddir}/src/libmps/libmps.la

%_moc.cpp: @srcdir@/%.h
	$(AM_V_GEN)
	QT_SELECT=${QT_VERSION} $(AM_V_at) $(MOC) $< -o $@

ui_%.h: @srcdir@/%.ui
	$(AM_V_GEN) 
	QT_SELECT=${QT_VERSION} $(AM_V_at) $(UIC) $< -o $@

%_qrc.cpp: @srcdir@/%.qrc
	$(AM_V_GEN)
	QT_SELECT=${QT_VERSION} $(AM_V_at) $(RCC) $< -o $@

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
	polynomialsolver.cpp  \
	rootsmodel.cpp  \
	polsyntaxhighlighter.cpp  \
	mpsolveworker.cpp

desktop_pure_source_files = \
	mainwindow.cpp \
	polfileeditorwindow.cpp \
	polfileeditor.cpp \
	qrootsrenderer.cpp

qml_pure_source_files = \
	qquickrootsrenderer.cpp \
	mainqmlview.cpp 

other_source_files = \
	main.cpp \
	rootsrenderer.cpp \
	root.cpp 

ui_files = \
	ui_mainwindow.h \
	ui_polfileeditor.h \
	ui_polfileeditorwindow.h 

ui_sfiles = \
	mainwindow.ui \
	polfileeditor.ui \
	polfileeditorwindow.ui

qrc_files = \
	resources.qrc

headers = \
	monomial.h \
	rootsrenderer.h  \
	rootsmodel.h  \
	root.h  \
	polsyntaxhighlighter.h  \
	polynomialsolver.h  \
	mpsolveworker.h

desktop_headers = \
	polfileeditor.h \
	polfileeditorwindow.h \
	mainwindow.h \
	qrootsrenderer.h

qml_headers = \
	mainqmlview.h \
	qquickrootsrenderer.h 

if BUILD_QML_UI
pure_source_files += $(qml_pure_source_files)
headers += $(qml_headers)
else
pure_source_files += $(desktop_pure_source_files)
headers += $(desktop_headers)
endif

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
	$(qml_headers) \
	$(desktop_headers) \
	xmpsolve.1 \
	xmpsolve.ico \
	xmpsolve.pro \
	icons/document-close.png \
	icons/document-open.png \
	icons/document-save.png \
	icons/document-new.png \
	icons/system-run.png \
	android/AndroidManifest.xml \
	android/res/drawable-mdpi/icon.png \
	android/res/drawable-hdpi/icon.png \
	android/res/drawable-ldpi/icon.png \
	Main.qml \
	LoadingIndicator.qml \
	MainView.qml \
	ApproximationList.qml \
	PolyInputField.qml \
	CombinedApproximationView.qml \
	SwitchableApproximationView.qml \
	loading.png 


BUILT_SOURCES = $(moc_files) \
	$(qrc_gen_files) 

if ! BUILD_QML_UI
BUILT_SOURCES += $(ui_files)
endif

endif # HAVE_UI
