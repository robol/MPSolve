#include "polfileeditordialog.h"
#include "ui_polfileeditordialog.h"
#include <QDebug>

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

void
PolFileEditorDialog::loadPolFile(QString polFilePath)
{
    m_polFilePath = polFilePath;

    QFile polFile(polFilePath);
    QString polFileContent;
    if (polFile.open(QIODevice::ReadOnly)) {
        QTextStream in(&polFile);
        polFileContent = in.readAll();
    }

    ui->polTextEdit->setPlainText(polFileContent);
    setWindowTitle(tr("Editing %1").arg(polFilePath));
}

void
PolFileEditorDialog::savePolFile()
{
    QFile outFile(m_polFilePath);
    if (outFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&outFile);
        out << ui->polTextEdit->toPlainText();
    }
}
