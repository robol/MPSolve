#ifndef POLFILEEDITORWINDOW_H
#define POLFILEEDITORWINDOW_H

#include <QMainWindow>
#include <QList>
#include "polfileeditor.h"

namespace Ui {
class PolFileEditorWindow;
}

namespace xmpsolve {

class PolFileEditorWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit PolFileEditorWindow(QWidget *parent = 0);
    ~PolFileEditorWindow();

    /**
     * @brief loadPolFile loads the file specified by path or simply
     * focus the tab containing it if it's already loaded.
     *
     * @brief path is the path to the file
     */
    void loadPolFile(QString path = QString());

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
     * @brief Close the given editor.
     */
    void closeEditor(PolFileEditor* editor);

    /**
     * @brief currentPolFile returns the path to the currently focused
     * .pol file.
     * @return A QString with the absolute path of the currently focused
     * .pol file.
     */
    QString currentPolFile();

    /**
     * @brief currentEditor can be used to access the current PolFileEditor focused
     * in the tabWidget.
     * @return A pointer to the focused instance of PolFileEditor, or NULL if there is none.
     */
    PolFileEditor* currentEditor();

signals:
    /**
     * @brief solvePoly is emitted when the user asks to solve a .pol file.
     * @param path is the path to the content of the file that the user wants to solve.
     */
    void solvePoly(QString content);

public slots:
    /**
     * @brief onEditorFilenameChanged handle the change of filename inside
     * and editor tab.
     */
    void onEditorFilenameChanged(QString);

    /**
     * @brief onEditorStateChanged handle the state changed of the editor tab.
     */
    void onEditorStateChanged(PolFileEditor::State);
    
private slots:
    void on_actionOpen_pol_file_triggered();

    void on_actionSave_triggered();

    void on_actionSolve_triggered();

    void on_actionClose_triggered();

    void on_actionClose_editor_triggered();

    void on_actionNew_triggered();

private:
    Ui::PolFileEditorWindow *ui;

    /**
     * @brief m_polFileEditors is a HashTable containing the associations
     * between paths on the filesystem and tab of the tabWidget.
     *
     * This allows to not open more than one tab at a time for a given file.
     */
    QMap<QString, PolFileEditor*> m_polFileEditors;

    /**
     * @brief Overriden closeEvent that handle the saving of files that
     * haven't been closed.
     */
    void closeEvent(QCloseEvent *);

    /**
     * @brief Close all opened tabs and prompt the user to save the unsaved ones.
     */
    void closeOpenedTabs();

    void showEvent(QShowEvent *event);

    /**
     * @brief setupIcons loads the icons needed by the Window.
     */
    void setupIcons();
};

} // End of namespace xmpsolve

#endif // POLFILEEDITORWINDOW_H
