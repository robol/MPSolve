#include "polsyntaxhighlighter.h"
#include <QRegExp>
#include <QTextCharFormat>
#include <QDebug>

PolSyntaxHighlighter::PolSyntaxHighlighter(QObject *parent) :
    QSyntaxHighlighter(parent)
{
}

void
PolSyntaxHighlighter::highlightBlock(const QString &text)
{
    QTextCharFormat commentFormat;
    commentFormat.setForeground(Qt::darkRed);

    QTextCharFormat statementFormat;
    statementFormat.setForeground(Qt::darkBlue);

    QRegExp commentRegexp("!.*");
    int index = text.indexOf(commentRegexp);
    while (index >= 0) {
        int length = commentRegexp.matchedLength();
        setFormat(index, length, commentFormat);
        index = text.indexOf(commentRegexp, index + length);
    }

    QRegExp assignmentRegExp("[^=]+\\s*=\\s*[^;]+\\s*;");
    index = text.indexOf(assignmentRegExp);
    while (index >= 0) {
        int length = assignmentRegExp.matchedLength();
        setFormat(index, length, statementFormat);
        index = text.indexOf(assignmentRegExp, index + length);
    }

    QRegExp statementRegexp(".+\\s*;");
    index = text.indexOf(statementRegexp);
    while (index >= 0) {
        int length = statementRegexp.matchedLength();
        setFormat(index, length, statementFormat);
        index = text.indexOf(statementRegexp, index + length);
    }
}
