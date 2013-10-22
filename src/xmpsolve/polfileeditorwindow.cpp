#include <QDebug>
#include <QFileInfo>
#include <QFileDialog>
#include <QMessageBox>
#include "polfileeditorwindow.h"
#include "ui_polfileeditorwindow.h"
#include "mainwindow.h"

using namespace xmpsolve;

PolFileEditorWindow::PolFileEditorWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::PolFileEditorWindow)
{
    ui->setupUi(this);
}

void
PolFileEditorWindow::showEvent(QShowEvent *event)
{
    if (ui->tabWidget->count() == 0)
        loadPolFile("");
}

void
PolFileEditorWindow::loadPolFile(QString path)
{
    PolFileEditor * editor = NULL;
    int currentIndex = 0;

    // We have 2 possibilities here:
    // - path is a path to a file, and we open it
    // - path is empty, and we have to create a new file
    if (! path.isEmpty()) {
        // Get the absolute path of the file to ensure an unique representation of each one.
        QFileInfo info(path);
        path = info.absoluteFilePath();
        editor =  m_polFileEditors.value(path, NULL);
    }

    if (editor == NULL) {
        editor = new PolFileEditor(this, path);

        QObject::connect(editor, SIGNAL(filenameChanged(QString)),
                         this, SLOT(onEditorFilenameChanged(QString)));

        if (! path.isEmpty())
            m_polFileEditors.insert(path, editor);

        ui->tabWidget->insertTab(0, editor, QIcon::fromTheme("text"),
                                 path.isEmpty() ? "New file" : path.split("/").last());
        QObject::connect(editor, SIGNAL(stateChanged(PolFileEditor::State)),
                         this, SLOT(onEditorStateChanged(PolFileEditor::State)));
    }
    else {
        currentIndex = ui->tabWidget->indexOf(editor);
    }

    ui->tabWidget->setCurrentIndex(currentIndex);
}

void
PolFileEditorWindow::savePolFile()
{
    PolFileEditor * editor = currentEditor();
    editor->savePolFile();
}

void
PolFileEditorWindow::closePolFile(QString path)
{
    PolFileEditor * editor = m_polFileEditors.value(path, NULL);
    QString filename = path.split("/").last();

    if (editor != NULL) {

        if (editor->state() == PolFileEditor::MODIFIED) {
            QMessageBox::StandardButton reply = QMessageBox::question(this,
                tr("Save changes to %1 before closing?").arg(filename),
                tr("%1 has been edited but no saved. Do you want to save it before closing it?").arg(
                    filename),
                QMessageBox::Yes|QMessageBox::No);

            if (reply == QMessageBox::Yes) {
                editor->savePolFile();
            }
        }

        int editorIndex = ui->tabWidget->indexOf(editor);
        ui->tabWidget->removeTab(editorIndex);
        m_polFileEditors.remove(path);

        delete editor;
    }

    if (ui->tabWidget->count() == 0) {
        loadPolFile("");
    }
}

PolFileEditorWindow::~PolFileEditorWindow()
{
    // Fake a close button press.
    on_actionClose_editor_triggered();

    delete ui;
}

void PolFileEditorWindow::onEditorFilenameChanged(QString newName)
{
    PolFileEditor *editor = static_cast<PolFileEditor*> (sender());
    int editorIndex = ui->tabWidget->indexOf(editor);

    ui->tabWidget->setTabText(editorIndex, newName.split("/").last());

    // Find the old key and remove it. Then insert the new one.
    QString oldKey = m_polFileEditors.key(editor, "");
    if (! oldKey.isEmpty()) {
        m_polFileEditors.remove(oldKey);
    }

    m_polFileEditors.insert(newName, editor);
}

void PolFileEditorWindow::onEditorStateChanged(PolFileEditor::State state)
{
    PolFileEditor *editor = static_cast<PolFileEditor*> (sender());
    QString path = editor->currentPolFile();

    if (path.isEmpty()) {
        path = "New file";
    }
    else {
        path = path.split("/").last();
    }

    int currentIndex = ui->tabWidget->indexOf(editor);

    switch (state) {
        case PolFileEditor::SAVED:
            ui->tabWidget->setTabText(currentIndex, path);
            break;
         case PolFileEditor::MODIFIED:
            ui->tabWidget->setTabText(currentIndex, path + "*");
            break;
    }
}

PolFileEditor* PolFileEditorWindow::currentEditor()
{
    PolFileEditor * editor = static_cast<PolFileEditor*> (ui->tabWidget->currentWidget());
    return editor;
}

QString PolFileEditorWindow::currentPolFile()
{
    PolFileEditor * editor = currentEditor();
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
    solvePoly(currentEditor()->content());
}

void PolFileEditorWindow::on_actionClose_triggered()
{
    closePolFile(currentPolFile());
}

void PolFileEditorWindow::on_actionNew_triggered()
{
    ui->tabWidget->insertTab(0, new PolFileEditor(this), QIcon::fromTheme("text"), "New file");
}

void PolFileEditorWindow::closeOpenedTabs()
{
    // Close the pol files before so we can get notifications about
    // unsaved files and so.
    foreach (PolFileEditor* editor, m_polFileEditors) {
        closePolFile(editor->currentPolFile());
    }

    // We need to handle the opened tabs that have not been
    // saved yet.
    while (ui->tabWidget->count() > 0) {

        ui->tabWidget->setCurrentIndex(0);
        PolFileEditor *editor = static_cast<PolFileEditor*>(ui->tabWidget->widget(0));

        if (editor->state() == PolFileEditor::MODIFIED) {

            int reply = QMessageBox::question(this, tr("Save changes to this file?"),
                                              tr("This file has not been saved. Do you want to do it now?"),
                                              QMessageBox::No, QMessageBox::Yes);
            if (reply == QMessageBox::Yes) {
                QString path = QFileDialog::getSaveFileName(this, "Save .pol file", "", "Pol files (*.pol)");
                if (! path.isEmpty()) {
                    editor->savePolFile(path);
                }
            }

        }

        ui->tabWidget->removeTab(0);
        delete editor;
    }
}

void PolFileEditorWindow::on_actionClose_editor_triggered()
{
    closeOpenedTabs();
    close();
}

void PolFileEditorWindow::closeEvent(QCloseEvent *event)
{
    closeOpenedTabs();
    close();
}
