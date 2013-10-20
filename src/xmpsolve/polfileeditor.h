#ifndef POLFILEEDITOR_H
#define POLFILEEDITOR_H

#include <QWidget>
#include "polsyntaxhighlighter.h"

namespace Ui {
class PolFileEditor;
}

class PolFileEditor : public QWidget
{
    Q_OBJECT
    
public:
    explicit PolFileEditor(QWidget *parent = 0, QString path = QString());
    void savePolFile();

    /**
     * @brief currentPolFile returns the currently opened .pol file
     * @return a QString with the absolute path of the file
     */
    QString currentPolFile();

    ~PolFileEditor();
    
private:
    Ui::PolFileEditor *ui;
    PolSyntaxHighlighter *m_syntaxHighlighter;
    QString m_polFilePath;

};

#endif // POLFILEEDITOR_H
