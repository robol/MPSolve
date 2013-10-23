#include "polfileeditor.h"
#include "ui_polfileeditor.h"
#include <QTextStream>
#include <QFileInfo>
#include <QFileDialog>
#include <QDebug>

using namespace xmpsolve;

PolFileEditor::PolFileEditor(QWidget *parent, QString path) :
    QWidget(parent),
    ui(new Ui::PolFileEditor)
{
    ui->setupUi(this);

    m_polFilePath = path;

    QString polFileContent;

    if (! path.isEmpty()) {
        QFile polFile(m_polFilePath);

        if (polFile.open(QIODevice::ReadOnly)) {
            QTextStream in(&polFile);
            polFileContent = in.readAll();
        }
    }

    m_syntaxHighlighter = new PolSyntaxHighlighter(ui->plainTextEdit->document());
    ui->plainTextEdit->setPlainText(polFileContent);
    m_syntaxHighlighter->setDocument(ui->plainTextEdit->document());

    QObject::connect(ui->plainTextEdit, SIGNAL(modificationChanged(bool)),
                     this, SLOT(onTextEditChanged(bool)));
}

void PolFileEditor::savePolFile(QString path)
{
    if (! path.isEmpty ())
        m_polFilePath = path;

    if (m_polFilePath.isEmpty()) {
        m_polFilePath = QFileDialog::getSaveFileName(this, "Save .pol file", "", "Pol files (*.pol)");
        filenameChanged(m_polFilePath);
    }

    QFile outFile(m_polFilePath);

    if (outFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&outFile);
        out << ui->plainTextEdit->toPlainText();
    }

    m_state = SAVED;
    ui->plainTextEdit->document()->setModified(false);
    stateChanged(m_state);
}

bool PolFileEditor::isEmpty()
{
    return ui->plainTextEdit->toPlainText() == "";
}

QString PolFileEditor::content()
{
    return ui->plainTextEdit->toPlainText();
}

PolFileEditor::State PolFileEditor::state()
{
    return m_state;
}

void PolFileEditor::onTextEditChanged(bool modified)
{
    m_state = modified ? MODIFIED : SAVED;
    stateChanged(m_state);
}

QString PolFileEditor::currentPolFile()
{
    if (m_polFilePath.isEmpty())
        return m_polFilePath;
    else {
        QFileInfo info(m_polFilePath);
        return info.absoluteFilePath();
    }
}

PolFileEditor::~PolFileEditor()
{
    delete ui;
    delete m_syntaxHighlighter;
}
