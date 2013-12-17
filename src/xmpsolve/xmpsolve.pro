# Project file for xmpsolve. Please note that this file assumes that you have
# compiled MPSolve in mpsolve-x.y.z/_build. 
#
# To accomplish this you can just enter the mpsolve-x.y.z folder and then
# $ mkdir _build
# $ cd _build 
# $ ../configure && make

TEMPLATE = app

CONFIG += mps_qml
# CONFIG += mps_widgets

QT       += core gui
mps_widgets {
    greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
}
mps_qml {
    greaterThan(QT_MAJOR_VERSION, 4) QT += quick
}

#
# This section addresses Android specific configuration. This is a bit atypical,
# w.r.t other platforms, since you will probably have to compile a copy of GMP
# and a separate version of libmps. Those will need to be installed in a different
# directory and then xmpsolve will be linked against them.
#
android {

    CONFIG += mobility
    MOBILITY =

    # App configuration. TODO: It should be possible (in the future, not right now)
    # to specify an icon here.
    ANDROID_PACKAGE = it.unipi.dm.mpsolve
    ANDROID_APP_NAME = MPSolve

    # Customize this to match your current setup. This way the setup points to a directory inside
    # mpsolve-x.y.z. This setup can be obtained with the script in tools/android-build-libmps.sh
    ANDROID_ROOT = $${PWD}/../../android-ext

    # We need -DMPS_USE_BUILTIN_COMPLEX since Android uses tiny bionic without complex
    # arithmetic support.
    QMAKE_CXXFLAGS += -I$${ANDROID_ROOT}/include \
        -include $${PWD}/../../android-ext/tarballs/mpsolve-build/config.h

    # Link against locally compiled libmps and libgmp.
    LIBS += $${ANDROID_ROOT}/lib/libmps.a $${ANDROID_ROOT}/lib/libgmp.a

}

!android {
    QMAKE_CXXFLAGS += -include $${PWD}/../../_build/config.h
    INCLUDEPATH += $${PWD}/../../include/
    LIBS += $${PWD}/../../_build/src/libmps/.libs/libmps.so -lgmp
}

# Input
HEADERS += ./mpsolveworker.h \
           ./polynomialsolver.h \
           ./root.h \
           ./rootsrenderer.h \
	   ./polynomialparser.h \
	   ./monomial.h \
	   ./polynomial.h \
           ./rootsmodel.h \
           ./polsyntaxhighlighter.h

SOURCES += ./main.cpp \
           ./mpsolveworker.cpp \
           ./polynomialsolver.cpp \
           ./root.cpp \
           ./rootsrenderer.cpp \
           ./polynomialparser.cpp \
           ./monomial.cpp \
           ./polynomial.cpp \
           ./rootsmodel.cpp \
           ./polsyntaxhighlighter.cpp

RESOURCES += \
    resources.qrc

OTHER_FILES += \
    android/AndroidManifest.xml \
    Main.qml \
    ApproximationList.qml \
    PolyInputField.qml \
    SwitchableApproximationView.qml \
    CombinedApproximationView.qml \
    MainView.qml \
    LoadingIndicator.qml

ANDROID_PACKAGE_SOURCE_DIR = $$PWD/android

# QtWidgets configuration
mps_widgets {
    HEADERS += mainwindow.h \
               polfileeditor.h \
               polfileeditorwindow.h \
               qrootsrenderer.h

    FORMS += ./mainwindow.ui \
             ./polfileeditor.ui \
             ./polfileeditorwindow.ui

    SOURCES += mainwindow.cpp \
           ./polfileeditor.cpp \
           ./polfileeditorwindow.cpp \
           ./qrootsrenderer.cpp

}

# QtQuick configuration
mps_qml {
    QMAKE_CXXFLAGS += -DMPS_USE_QML

    HEADERS += qquickrootsrenderer.h  \
               mainqmlview.h

    SOURCES += ./mainqmlview.cpp \
           ./qquickrootsrenderer.cpp
}
