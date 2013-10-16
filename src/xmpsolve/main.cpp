#include <qglobal.h>

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
#include <QtGui/QApplication>
#else
#include <QtWidgets/QApplication>
#endif
#include "mainwindow.h"

using namespace xmpsolve;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainWindow w;
    w.show();

    // In case the user wants to open a .pol file, try load it
    if (argc > 1) {
        w.openPolFile(argv[1]);
    }
    
    return a.exec();
}
