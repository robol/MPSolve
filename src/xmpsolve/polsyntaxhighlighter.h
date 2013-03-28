#ifndef POLSYNTAXHIGHLIGHTER_H
#define POLSYNTAXHIGHLIGHTER_H

#include <QSyntaxHighlighter>

class PolSyntaxHighlighter : public QSyntaxHighlighter
{
    Q_OBJECT
public:
    explicit PolSyntaxHighlighter(QObject *parent = 0);

    void highlightBlock(const QString &text);
    
signals:
    
public slots:
    
};

#endif // POLSYNTAXHIGHLIGHTER_H
