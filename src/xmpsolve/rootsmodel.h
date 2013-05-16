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
        SHORT_APPROXIMATION,
        ROOT
    };

    explicit RootsModel(QObject *parent = 0);
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    QHash<int, QByteArray> roleNames() const;

    void setRoots(QList<Root*> roots);

private:
    QList<Root*> m_roots;
    int length;
    
signals:
    
public slots:
    
};

} // namespace xmpsolve

#endif // XMPSOLVE_ROOTSMODEL_H
