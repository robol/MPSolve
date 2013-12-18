#include <qglobal.h>

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
#include <QtGui/QApplication>
#else
#ifdef MPS_USE_QML
#include <QtGui/QGuiApplication>
#else
#include <QtWidgets/QApplication>
#endif
#endif

#ifdef MPS_USE_QML
#include <QtQml/QQmlApplicationEngine>
#include "mainqmlview.h"
#else
#include "mainwindow.h"
#endif

using namespace xmpsolve;

int main(int argc, char *argv[])
{
#ifdef MPS_USE_QML
    QGuiApplication a(argc, argv);
#else
    QApplication a(argc, argv);
#endif

#ifdef MPS_USE_QML
    MainQmlView w;
#else
    MainWindow w;

    // In case the user wants to open a .pol file, try to load it
    if (argc > 1) {
        w.openPolFile(argv[1]);
    }

    w.show();
#endif
    
    return a.exec();
}
