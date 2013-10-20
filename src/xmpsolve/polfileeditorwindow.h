#ifndef POLFILEEDITORWINDOW_H
#define POLFILEEDITORWINDOW_H

#include <QMainWindow>
#include <QList>
#include "polfileeditor.h"

namespace Ui {
class PolFileEditorWindow;
}

class PolFileEditorWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit PolFileEditorWindow(QWidget *parent = 0);
    ~PolFileEditorWindow();

    /**
     * @brief loadPolFile loads the file specified by path or simply
     * focus the tab containing it if it's already loaded.
     */
    void loadPolFile(QString path);

    /**
     * @brief savePolFile save the pol file in the currently selected tab.
     */
    void savePolFile();

    /**
     * @brief closePolFile closes the tab of the given .pol file
     * @param path is the absolute path to the file to close.
     */
    void closePolFile(QString path);

    /**
     * @brief currentPolFile returns the path to the currently focused
     * .pol file.
     * @return A QString with the absolute path of the currently focused
     * .pol file.
     */
    QString currentPolFile();

signals:
    /**
     * @brief solvePoly is emitted when the user asks to solve a .pol file.
     * @param path is the path to the file that the user wants to solve.
     */
    void solvePoly(QString path);
    
private slots:
    void on_actionOpen_pol_file_triggered();

    void on_actionSave_triggered();

    void on_actionSolve_triggered();

    void on_actionClose_triggered();

    void on_actionClose_editor_triggered();

private:
    Ui::PolFileEditorWindow *ui;

    /**
     * @brief m_polFileEditors is a HashTable containing the associations
     * between paths on the filesystem and tab of the tabWidget.
     *
     * This allows to not open more than one tab at a time for a given file.
     */
    QMap<QString, PolFileEditor*> m_polFileEditors;
};

#endif // POLFILEEDITORWINDOW_H
