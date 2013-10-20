#ifndef POLSYNTAXHIGHLIGHTER_H
#define POLSYNTAXHIGHLIGHTER_H

#include <QSyntaxHighlighter>

namespace xmpsolve {

class PolSyntaxHighlighter : public QSyntaxHighlighter
{
    Q_OBJECT
public:
    explicit PolSyntaxHighlighter(QObject *parent = 0);

    void highlightBlock(const QString &text);
    
signals:
    
public slots:
    
};

} // End of namespace xmpsolve

#endif // POLSYNTAXHIGHLIGHTER_H
