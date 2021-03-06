#ifndef XMPSOLVE_ROOTSMODEL_H
#define XMPSOLVE_ROOTSMODEL_H

#include <QAbstractListModel>
#include "root.h"

namespace xmpsolve {

class RootsModel : public QAbstractListModel
{
    Q_OBJECT


public:

    enum Roles {
        RADIUS = Qt::UserRole + 1,
        STATUS,
        SHORT_APPROXIMATION,
        ROOT,
        MARKED
    };

    explicit RootsModel(QObject *parent = 0);

    Q_INVOKABLE int rowCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;

    QHash<int, QByteArray> roleNames() const;

    void setRoots(QList<Root*> roots);

    /**
     * @brief markRoot can be used to highlight an approximation.
     * @param i is the root to highlight, or -1 to clear any previous
     * highlighting.
     */
    Q_INVOKABLE void markRoot(int i = -1);

    double getPointX(int i) { return m_roots[i]->get_real_part(); }
    double getPointY(int i) { return m_roots[i]->get_imag_part(); }

private:
    QList<Root*> m_roots;
    int m_length;
    int m_marked_root;
    
signals:
    
public slots:
    
};

} // namespace xmpsolve

#endif // XMPSOLVE_ROOTSMODEL_H
