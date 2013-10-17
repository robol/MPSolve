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

    QObject::connect(ui->listRootsView->selectionModel(), SIGNAL(selectionChanged(QItemSelection,QItemSelection)),
                     this, SLOT(onlistRootsView_selectionChanged(QItemSelection,QItemSelection)));
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

    PolynomialBasis basis = MONOMIAL;
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
    openPolFile(selectedFile);
}

void xmpsolve::MainWindow::openPolFile(QString path)
{
    QFile polFile(path);

    // Clean previous polynomials, if any
    ui->polyLineEdit->clear();

    // Select the polynomial
    if (polFile.exists()) {
        m_selectedPolFile = path;
        ui->selectedFileLabel->setText(path.split("/").last());

        ui->selectedFileLabel->setEnabled(true);
        ui->editPolFileButton->setEnabled(true);
    }

    // Switch the view to the correct tab
    ui->tabWidget->setCurrentIndex(1);
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

void xmpsolve::MainWindow::onlistRootsView_selectionChanged(QItemSelection,QItemSelection)
{
    QModelIndexList indexes = ui->listRootsView->selectionModel()->selection().indexes();

    switch (indexes.length())
    {
        case 0:
            m_solver.rootsModel()->markRoot();
            break;
        case 1:
            m_solver.rootsModel()->markRoot(indexes.at(0).row());
            break;
        default:
            qDebug() << "More than 1 index selected";
            break;
    }
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

void xmpsolve::MainWindow::on_actionOpen_pol_file_triggered()
{
    on_openPolFileButton_clicked();
}



void xmpsolve::MainWindow::on_actionQuit_triggered()
{
    QApplication::exit();
}

void xmpsolve::MainWindow::on_actionAbout_MPSolve_triggered()
{
    QMessageBox::about(this, tr("About ") + PACKAGE_STRING,
                       QString("<h1>%1</h1>").arg(PACKAGE_STRING) +
                       "MPSolve is free software released under the GNU General Public License 3.<br>" +
                       "Further documentation and a bug tracker are available at our " +
                       "<a href=\"http://mpsolve.dm.unipi.it/mpsolve/\">website</a>. " +
                       "<h3>Authors:</h3>" +
                       " - Dario A. Bini &lt;<a href=\"mailto:bini@dm.unipi.it\">bini@dm.unipi.it</a>&gt;<br>" +
                       " - Giuseppe Fiorentino &lt;<a href=\"mailto:fiorent@dm.unipi.it\">fiorent@dm.unipi.it</a>&gt;<br>" +
                       " - Leonardo Robol &lt;<a href=\"mailto:leonardo.robol@sns.it\">leonardo.robol@sns.it</a>&gt; <br>");
}
