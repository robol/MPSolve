#include <QtGui/QApplication>
#include "mainwindow.h"

using namespace xmpsolve;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainWindow w;
    w.show();
    
    return a.exec();
}
