#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <polynomialsolver.h>
#include <QDebug>
#include <QPainter>
#include <QMessageBox>
#include <QFileDialog>
#include "polfileeditordialog.h"

using namespace xmpsolve;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Connect the signal of computation finished to its callback in
    // here.
    QObject::connect(&m_solver, SIGNAL(solved()),
                     this, SLOT(polynomial_solved()));

    ui->listRootsView->setModel(m_solver.rootsModel());
    ui->graphicsView->setModel(m_solver.rootsModel());
}

MainWindow::~MainWindow()
{
    delete ui;
}

void
MainWindow::lockInterface()
{
    ui->solveButton->setEnabled(false);
    ui->polFileSolveButton->setEnabled(false);
}

void
MainWindow::unlockInterface()
{
    ui->solveButton->setEnabled(true);
    ui->polFileSolveButton->setEnabled(true);
}

void MainWindow::on_solveButton_clicked()
{
    lockInterface();

    mps_algorithm selected_algorithm = (ui->algorithmComboBox->currentIndex() == 0) ?
                MPS_ALGORITHM_SECULAR_GA : MPS_ALGORITHM_STANDARD_MPSOLVE;

    PolynomialBasis basis;
    switch (ui->basisComboBox->currentIndex()) {
        case 0:
            basis = MONOMIAL;
            break;
         case 1:
            basis = CHEBYSHEV;
            break;
    }

    ui->statusBar->showMessage(tr("Solving polynomial..."));
    if (m_solver.solvePoly(ui->polyLineEdit->toPlainText(), basis,
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
MainWindow::polynomial_solved()
{
    ui->statusBar->showMessage(tr("Polynomial solved in %1ms").arg(m_solver.CPUTime()));

    unlockInterface();
}

void xmpsolve::MainWindow::on_openPolFileButton_clicked()
{
    // If the user click on open polFile, we need to check that
    // he has selected a valid .pol file and - if that's the
    // case, solve the associated polynomial.
    QString selectedFile = QFileDialog::getOpenFileName(this,
                                                        tr("Select .pol file"),
                                                        QString(),
                                                        "Pol files (*.pol);;Text files (*.txt)");
    if (! selectedFile.isEmpty())
    {
        QFile polFile(selectedFile);

        // Clean previous polynomials, if any
        ui->polyLineEdit->clear();

        // Select the polynomial
        if (polFile.exists()) {
            m_selectedPolFile = selectedFile;
            ui->selectedFileLabel->setText(selectedFile.split("/").last());

            ui->selectedFileLabel->setEnabled(true);
            ui->editPolFileButton->setEnabled(true);
        }
    }
}

void xmpsolve::MainWindow::on_listRootsView_clicked(const QModelIndex &index)
{
    ui->approximationDetailLabel->setText(
                m_solver.rootsModel()->data(index, RootsModel::SHORT_APPROXIMATION).toString());
    ui->approximationRadiusLabel->setText(
                m_solver.rootsModel()->data(index, RootsModel::RADIUS).toString());
    ui->approximationStatusLabel->setText(
                m_solver.rootsModel()->data(index, RootsModel::STATUS).toString());
}

void xmpsolve::MainWindow::on_editPolFileButton_clicked()
{
    PolFileEditorDialog dialog; // = new PolFileEditorDialog();
    dialog.loadPolFile(m_selectedPolFile);

    // Check if need to save the file.
    if (dialog.exec() == QDialog::Accepted) {
        dialog.savePolFile();
    }
}

void xmpsolve::MainWindow::on_polFileSolveButton_clicked()
{
    lockInterface();

    mps_algorithm selected_algorithm = (ui->algorithmComboBox->currentIndex() == 0) ?
                MPS_ALGORITHM_SECULAR_GA : MPS_ALGORITHM_STANDARD_MPSOLVE;

    // If we have a .pol file active, solve that. In the other case, try
    // to parse the input written by the user.
    if (! m_selectedPolFile.isEmpty()) {
        ui->statusBar->showMessage(tr("Solving %1...").arg(m_selectedPolFile));

        if (m_solver.solvePolFile(m_selectedPolFile, selected_algorithm,
                                  ui->digitsSpinBox->value()) < 0) {
            ui->statusBar->showMessage(tr("Polynomial parsing failed"));
            QMessageBox mbox(QMessageBox::Critical, tr("Error while parsing the polynomial"),
                             tr("The parser reported the following error: ") +
                             m_solver.errorMessage(), QMessageBox::Ok);
            mbox.exec();
            qDebug() << m_solver.errorMessage();
            unlockInterface();
        }
    }
    else {
        ui->statusBar->showMessage(tr("Please select a .pol file"));
        QMessageBox mbox(QMessageBox::Critical, tr("No .pol file selected"),
                         tr("Please select a valid .pol file to solve."),
                         QMessageBox::Ok);
        mbox.exec();
        unlockInterface();
    }
}
