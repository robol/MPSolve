#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <polynomialsolver.h>
#include <QDebug>
#include <QPainter>
#include <QMessageBox>

using namespace xmpsolve;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Connect the signal of computation finished to its callback in
    // here.
    QObject::connect(&m_solver, SIGNAL(solved(QList<Root*>)),
                     this, SLOT(polynomial_solved(QList<Root*>)));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void
MainWindow::lockInterface()
{
    ui->plainTextEdit->setEnabled(false);
    ui->solveButton->setEnabled(false);
}

void
MainWindow::unlockInterface()
{
    ui->plainTextEdit->setEnabled(true);
    ui->solveButton->setEnabled(true);
}

void MainWindow::on_solveButton_clicked()
{
    ui->statusBar->showMessage(tr("Solving polynomial..."));
    lockInterface();

    mps_algorithm selected_algorithm = (ui->algorithmComboBox->currentIndex() == 0) ?
                MPS_ALGORITHM_SECULAR_GA : MPS_ALGORITHM_STANDARD_MPSOLVE;

    if (m_solver.solvePoly(ui->polyLineEdit->toPlainText(),
                          selected_algorithm) < 0)
    {
        ui->statusBar->showMessage(tr("Polynomial parsing failed"));
        QMessageBox mbox(QMessageBox::Critical, tr("Error while parsing the polynomial"),
                         tr("The parser reported the following error: ") +
                         m_solver.errorMessage(), QMessageBox::Ok);
        mbox.exec();
        unlockInterface();
    }
}

void
MainWindow::polynomial_solved(QList<Root*> roots)
{
    QString output;
    ui->statusBar->showMessage(tr("Polynomial solved in %1ms").arg(m_solver.CPUTime()));

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

    unlockInterface();

    ui->plainTextEdit->clear();
    ui->plainTextEdit->insertPlainText(output);

    // Draw the result on the graphicsView.
    ui->graphicsView->setRoots(roots);
}
