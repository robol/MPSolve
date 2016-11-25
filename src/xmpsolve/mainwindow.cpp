#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <polynomialsolver.h>
#include <QDebug>
#include <QPainter>
#include <QMessageBox>
#include <QFileDialog>
#include "polfileeditorwindow.h"

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

    m_polFileEditorWindow = NULL;

    // Synchronize the selection of the roots with the rootsView, so we can
    // focus the currently selected root.
    QObject::connect(ui->listRootsView->selectionModel(),
		     SIGNAL(selectionChanged(QItemSelection,QItemSelection)),
                     this,
		     SLOT(onlistRootsView_selectionChanged(QItemSelection,QItemSelection)));
}

MainWindow::~MainWindow()
{
    if (m_polFileEditorWindow != NULL)
        delete m_polFileEditorWindow;
    delete ui;
}

void
MainWindow::openEditor(QString polFile)
{
    if (m_polFileEditorWindow == NULL) {
        m_polFileEditorWindow = new PolFileEditorWindow();

        // Solve .pol files requested by the user
        QObject::connect(m_polFileEditorWindow, SIGNAL(solvePoly(QString)),
                         this, SLOT(onSolvePolFileRequested(QString)));

        QObject::connect(m_polFileEditorWindow, SIGNAL(destroyed()),
                         this, SLOT(onPolFileEditorWindowDestroyed()));
    }

    if (! polFile.isEmpty()) {
        m_polFileEditorWindow->loadPolFile(polFile);
    }

    m_polFileEditorWindow->raise();
    m_polFileEditorWindow->show();
}

void
MainWindow::onPolFileEditorWindowDestroyed()
{
    m_polFileEditorWindow->deleteLater();
    m_polFileEditorWindow = NULL;
}

void
MainWindow::lockInterface()
{
    ui->solveButton->setEnabled(false);
}

void
MainWindow::unlockInterface()
{
    ui->solveButton->setEnabled(true);
}

int MainWindow::requiredDigits()
{
    return ui->digitsSpinBox->value();
}

mps_algorithm MainWindow::selectedAlgorithm()
{
    mps_algorithm selected_algorithm = (ui->algorithmComboBox->currentIndex() == 0) ?
                MPS_ALGORITHM_SECULAR_GA : MPS_ALGORITHM_STANDARD_MPSOLVE;
    return selected_algorithm;
}

PolynomialBasis MainWindow::polynomialBasis()
{
    PolynomialBasis basis = MONOMIAL;
    switch (ui->basisComboBox->currentIndex()) {
        case 0:
            basis = MONOMIAL;
            break;
         case 1:
            basis = CHEBYSHEV;
            break;
    }

    return basis;
}

mps_output_goal MainWindow::selectedGoal()
{
    switch (ui->goalComboBox->currentIndex()) {
    case 0:
        return MPS_OUTPUT_GOAL_ISOLATE;
        break;

    case 1:
        return MPS_OUTPUT_GOAL_APPROXIMATE;
        break;
    default:
        // This is kept for now as a default value, but this
        // point should not be reached.
        return MPS_OUTPUT_GOAL_APPROXIMATE;
    }
}

void MainWindow::on_solveButton_clicked()
{
    lockInterface();

    ui->statusBar->showMessage(tr("Solving polynomial..."));
    if (m_solver.solvePoly(ui->polyLineEdit->toPlainText(), polynomialBasis(),
                           selectedAlgorithm(), requiredDigits(),
                           selectedGoal()) < 0)
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

void xmpsolve::MainWindow::openPolFile(QString path)
{
    openEditor(path);
}

void xmpsolve::MainWindow::onSolvePolFileRequested(QString content)
{
    lockInterface();

    if (m_solver.solvePolFileFromContent(content, selectedAlgorithm(),
                                         requiredDigits(), selectedGoal()) < 0) {
        ui->statusBar->showMessage(tr("Polynomial parsing failed"));

        QMessageBox::critical(this,
                              tr("Error while parsing the polynomial"),
                              tr("The parser reported the following error: ") +
                              m_solver.errorMessage(), QMessageBox::Ok);

        unlockInterface();
    }

    // Grab focus so we can display the result.
    raise();
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

void xmpsolve::MainWindow::on_actionOpen_pol_file_triggered()
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

void xmpsolve::MainWindow::closeEvent(QCloseEvent *)
{
    QApplication::exit();
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

void xmpsolve::MainWindow::on_actionOpen_editor_triggered()
{
    openEditor();
}

void xmpsolve::MainWindow::on_actionAbort_computations_triggered()
{
    m_solver.abortComputations();
    ui->statusBar->showMessage(tr("Waiting for MPSolve to complete the current operation..."));
}

void xmpsolve::MainWindow::on_openPolFileButton_clicked()
{
    on_actionOpen_pol_file_triggered();
}
