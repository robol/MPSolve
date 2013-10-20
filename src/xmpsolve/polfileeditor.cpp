#include "polfileeditor.h"
#include "ui_polfileeditor.h"
#include <QTextStream>
#include <QFileInfo>

PolFileEditor::PolFileEditor(QWidget *parent, QString path) :
    QWidget(parent),
    ui(new Ui::PolFileEditor)
{
    ui->setupUi(this);

    m_polFilePath = path;
    QFile polFile(m_polFilePath);
    QString polFileContent;

    if (polFile.open(QIODevice::ReadOnly)) {
        QTextStream in(&polFile);
        polFileContent = in.readAll();
    }

    m_syntaxHighlighter = new PolSyntaxHighlighter(ui->plainTextEdit->document());
    ui->plainTextEdit->setPlainText(polFileContent);
    m_syntaxHighlighter->setDocument(ui->plainTextEdit->document());
}

void PolFileEditor::savePolFile()
{
    QFile outFile(m_polFilePath);

    if (outFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&outFile);
        out << ui->plainTextEdit->toPlainText();
    }
}

QString PolFileEditor::currentPolFile()
{
    QFileInfo info(m_polFilePath);
    return info.absoluteFilePath();
}

PolFileEditor::~PolFileEditor()
{
    delete ui;
    delete m_syntaxHighlighter;
}
