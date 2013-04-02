#include "polfileeditordialog.h"
#include "ui_polfileeditordialog.h"
#include "polsyntaxhighlighter.h"
#include <QDebug>
#include <QFile>

PolFileEditorDialog::PolFileEditorDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PolFileEditorDialog)
{
    ui->setupUi(this);
    m_highlighter = new PolSyntaxHighlighter(ui->polTextEdit->document());
}

PolFileEditorDialog::~PolFileEditorDialog()
{
    delete ui;
    delete m_highlighter;
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
    m_highlighter->setDocument(ui->polTextEdit->document());
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
