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

    /**
     * @brief State of the document
     */
    enum State {
        SAVED,
        MODIFIED
    };

    explicit PolFileEditor(QWidget *parent = 0, QString path = QString());
    void savePolFile();

    /**
     * @brief currentPolFile returns the currently opened .pol file
     * @return a QString with the absolute path of the file
     */
    QString currentPolFile();

    /**
     * @brief isEmpty returns true if the editor has no content in it.
     * @return
     */
    bool isEmpty();

    /**
     * @brief state returns the current state of the document
     * @return the value of m_state.
     */
    State state();

    /**
     * @brief content returns the content of the editor
     * @return A QString with the current content of the Editor.
     */
    QString content();

    ~PolFileEditor();

signals:
    void filenameChanged(QString filename);
    void stateChanged(PolFileEditor::State);

private slots:
    void onTextEditChanged(bool);
    
private:
    Ui::PolFileEditor *ui;
    PolSyntaxHighlighter *m_syntaxHighlighter;
    QString m_polFilePath;
    State m_state;

};

#endif // POLFILEEDITOR_H
