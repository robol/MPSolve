#include "polfileeditordialog.h"
#include "ui_polfileeditordialog.h"

PolFileEditorDialog::PolFileEditorDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PolFileEditorDialog)
{
    ui->setupUi(this);
}

PolFileEditorDialog::~PolFileEditorDialog()
{
    delete ui;
}
