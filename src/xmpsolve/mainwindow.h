#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QItemSelection>
#include "root.h"
#include "polynomialsolver.h"
#include "polfileeditorwindow.h"
#include <mps/mps.h>

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

    /**
     * @brief selectedAlgorithm return the algorithm currently selected by the user.
     * @return An mps_algorithm with the value selected by the user.
     */
    mps_algorithm selectedAlgorithm();

    /**
     * @brief requiredDigits returns the number of accuracy digits required by the user
     * @return A positive integer
     */
    int requiredDigits();

    /**
     * @brief selectedGoal returns the goal for the computation. Can be either set to
     * MPS_OUTPUT_GOAL_ISOLATE or to MPS_OUTPUT_GOAL_APPROXIMATE.
     *
     * @return The goal for the computation selected by the user.
     */
    mps_output_goal selectedGoal();

    /**
     * @brief polynomialBasis returns the polynomial basis selected by the user.
     * @return A PolynomialBasis value.
     */
    PolynomialBasis polynomialBasis();

    /**
     * @brief openEditor performs the necessary steps to set up the PolFileEditorWindow
     * and activate it.
     *
     * @param polFile If a non empty string is passed as argument, the file is opened
     * in the Editor.
     *
     * In the case the window is already available, it just focuses it.
     */
    void openEditor(QString polFile = "");

public slots:
    void polynomial_solved();

    /**
     * @brief openPolFile loads a .pol file given its path
     * @param path The path to the .pol file
     */
    void openPolFile(QString path);

    void onlistRootsView_selectionChanged(QItemSelection, QItemSelection);
    void onSolvePolFileRequested(QString path);
    
private slots:
    void on_solveButton_clicked();
    void lockInterface();
    void unlockInterface();

    void on_listRootsView_clicked(const QModelIndex &index);

    void on_actionOpen_pol_file_triggered();

    void on_actionQuit_triggered();

    void on_actionAbout_MPSolve_triggered();

    void on_actionOpen_editor_triggered();

    void on_actionAbort_computations_triggered();

    void closeEvent(QCloseEvent *);

    void onPolFileEditorWindowDestroyed();

    void on_openPolFileButton_clicked();

    void on_zoomInButton_clicked();

    void on_zoomOutButton_clicked();

private:
    Ui::MainWindow *ui;
    PolynomialSolver m_solver;
    QString m_selectedPolFile;
    PolFileEditorWindow *m_polFileEditorWindow;
};

} // Namespace xmpsolve

#endif // MAINWINDOW_H
