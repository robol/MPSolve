#include <QDebug>
#include <QFileInfo>
#include <QFileDialog>
#include "polfileeditorwindow.h"
#include "ui_polfileeditorwindow.h"
#include "mainwindow.h"

PolFileEditorWindow::PolFileEditorWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::PolFileEditorWindow)
{
    ui->setupUi(this);
}

void
PolFileEditorWindow::loadPolFile(QString path)
{
    // Get the absolute path of the file to ensure an unique representation
    // of each one.
    QFileInfo info(path);
    path = info.absoluteFilePath();

    PolFileEditor * editor = m_polFileEditors.value(path, NULL);
    int currentIndex = 0;

    if (editor == NULL) {
        editor = new PolFileEditor(this, path);
        m_polFileEditors.insert(path, editor);

        ui->tabWidget->insertTab(0, editor, QIcon::fromTheme("text"),
                                 path.split("/").last());
    }
    else {
        currentIndex = ui->tabWidget->indexOf(editor);
    }

    ui->tabWidget->setCurrentIndex(currentIndex);
}

void
PolFileEditorWindow::savePolFile()
{
    PolFileEditor * editor = static_cast<PolFileEditor*> (ui->tabWidget->currentWidget());
    editor->savePolFile();
}

void
PolFileEditorWindow::closePolFile(QString path)
{
    PolFileEditor * editor = m_polFileEditors.value(path, NULL);

    if (editor != NULL) {

        // TODO: We may want to handle the "You haven't saved the file..."
        // kind of notifications here.

        int editorIndex = ui->tabWidget->indexOf(editor);
        ui->tabWidget->removeTab(editorIndex);
        m_polFileEditors.remove(path);

        delete editor;
    }
}

PolFileEditorWindow::~PolFileEditorWindow()
{
    delete ui;
}

QString PolFileEditorWindow::currentPolFile()
{
    PolFileEditor * editor = static_cast<PolFileEditor*> (ui->tabWidget->currentWidget());
    if (editor != NULL)
        return editor->currentPolFile();
    else
        return QString();
}

void PolFileEditorWindow::on_actionOpen_pol_file_triggered()
{
    QString path = QFileDialog::getOpenFileName(this, "Select .pol file", "", "Pol files (*.pol)");
    loadPolFile(path);
}

void PolFileEditorWindow::on_actionSave_triggered()
{
    savePolFile();
}

void PolFileEditorWindow::on_actionSolve_triggered()
{
    solvePoly(currentPolFile());
}

void PolFileEditorWindow::on_actionClose_triggered()
{
    closePolFile(currentPolFile());
}

void PolFileEditorWindow::on_actionClose_editor_triggered()
{
    close();
}
