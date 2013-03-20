#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <polynomialsolver.h>
#include <QDebug>
#include <QPainter>
#include <QMessageBox>

#include <gmpxx.h>

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
                          selected_algorithm, ui->digitsSpinBox->value()) < 0)
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
        int digits = ui->digitsSpinBox->value();
        char * buffer = new char[2 * digits + 15];

        gmp_sprintf (buffer, "%.*Ff + %.*Ffi", digits, mpc_Re (roots[i]->value),
                     digits, mpc_Im (roots[i]->value));

        output.append(buffer);

        delete [] buffer;
        output.append("\n");
    }

    unlockInterface();

    ui->plainTextEdit->clear();
    ui->plainTextEdit->insertPlainText(output);

    // Draw the result on the graphicsView.
    ui->graphicsView->setRoots(roots);
}
