#ifndef POLFILEEDITORDIALOG_H
#define POLFILEEDITORDIALOG_H

#include <QDialog>
#include "polsyntaxhighlighter.h"

namespace Ui {
class PolFileEditorDialog;
}

class PolFileEditorDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit PolFileEditorDialog(QWidget *parent = 0);
    ~PolFileEditorDialog();

    void loadPolFile(QString polFilePath);
    void savePolFile();
    
private:
    Ui::PolFileEditorDialog *ui;
    PolSyntaxHighlighter *m_highlighter;
    QString m_polFilePath;
};

#endif // POLFILEEDITORDIALOG_H
