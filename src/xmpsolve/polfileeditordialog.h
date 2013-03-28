#ifndef POLFILEEDITORDIALOG_H
#define POLFILEEDITORDIALOG_H

#include <QDialog>

namespace Ui {
class PolFileEditorDialog;
}

class PolFileEditorDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit PolFileEditorDialog(QWidget *parent = 0);
    ~PolFileEditorDialog();
    
private:
    Ui::PolFileEditorDialog *ui;
};

#endif // POLFILEEDITORDIALOG_H
