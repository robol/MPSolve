#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "root.h"
#include "polynomialsolver.h"
#include<mps/mps.h>

namespace Ui {
    class MainWindow;
}

namespace xmpsolve {

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void polynomial_solved(QList<Root*>);
    
private slots:
    void on_solveButton_clicked();
    void lockInterface();
    void unlockInterface();

private:
    Ui::MainWindow *ui;
    PolynomialSolver m_solver;
};

} // Namespace xmpsolve

#endif // MAINWINDOW_H
