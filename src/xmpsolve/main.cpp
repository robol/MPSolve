#include <qglobal.h>

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
#include <QtGui/QApplication>
#else
#include <QtWidgets/QApplication>
#include <QtQuick/QQuickView>
#endif
#include "mainwindow.h"

#ifdef MPS_USE_QML
#include "mainqmlview.h"
#endif

using namespace xmpsolve;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

#ifdef MPS_USE_QML
    MainQmlView w;
#else
    MainWindow w;

    // In case the user wants to open a .pol file, try load it
    if (argc > 1) {
        w.openPolFile(argv[1]);
    }
#endif

    w.show();
    
    return a.exec();
}
