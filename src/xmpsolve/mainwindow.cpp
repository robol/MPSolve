#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <polynomialsolver.h>
#include <QDebug>
#include <QPainter>

using namespace xmpsolve;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_solveButton_clicked()
{
    /* Create our polynomial solver */
    PolynomialSolver *solver = new PolynomialSolver();

    // Connect the signal of computation finished to its callback in
    // here.
    QObject::connect(solver, SIGNAL(solved(QList<Root*>)),
                     this, SLOT(polynomial_solved(QList<Root*>)));


    ui->statusBar->showMessage(tr("Solving polynomial..."));

    if (solver->solvePoly(ui->polyLineEdit->text()) < 0)
    {
            QList<Root*> empty_list;
            polynomial_solved(empty_list);
    }
}

void
MainWindow::polynomial_solved(QList<Root*> roots)
{
    QString output;
    ui->statusBar->showMessage(tr("Polynomial solved"));

    for (int i = 0; i < roots.length(); i++)
    {
        // output.append (QString ("Root %1: ").arg(i));
        if (roots[i]->get_imag_part() != 0)
            output.append (QString ("%1 + %2i").arg(roots[i]->get_real_part()).arg(roots[i]->get_imag_part()));
        else
            output.append (QString ("%1").arg(roots[i]->get_real_part()));

        // output.append (QString("Radius: %1\n").arg(roots[i]->get_radius()));
        output.append("\n");
    }

    ui->plainTextEdit->clear();
    ui->plainTextEdit->insertPlainText(output);

    // Draw the result on the graphicsView.
    ui->graphicsView->setRoots(roots);
}
