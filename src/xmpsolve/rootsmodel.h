#ifndef XMPSOLVE_ROOTSMODEL_H
#define XMPSOLVE_ROOTSMODEL_H

#include <QAbstractListModel>
#include "root.h"

namespace xmpsolve {

class RootsModel : public QAbstractListModel
{
    Q_OBJECT


public:

    enum RootsModelRoles {
        RADIUS = Qt::UserRole + 1,
        STATUS,
        SHORT_APPROXIMATION
    };

    explicit RootsModel(QObject *parent = 0);
    int rowCount(const QModelIndex &parent) const;
    QVariant data(const QModelIndex &index, int role) const;

    void setRoots(QList<Root*> roots);

private:
    QList<Root*> m_roots;
    int length;
    
signals:
    
public slots:
    
};

} // namespace xmpsolve

#endif // XMPSOLVE_ROOTSMODEL_H
