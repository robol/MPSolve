#ifndef QQUICKROOTSRENDERER_H
#define QQUICKROOTSRENDERER_H

#include <QQuickPaintedItem>
#include "rootsmodel.h"
#include "rootsrenderer.h"

namespace xmpsolve {

    class QQuickRootsRenderer : public QQuickPaintedItem, public RootsRenderer
    {
        Q_OBJECT
        Q_PROPERTY(QVariant model READ model WRITE setModel)

    public:
        explicit QQuickRootsRenderer(QQuickItem *parent = 0);

        void paint(QPainter * painter);

        QVariant model() { return QVariant::fromValue(m_model); }void
        setModel(QVariant model);

    signals:

    public slots:
        void reloadRootsWrapper();

    };

}

#endif // QQUICKROOTSRENDERER_H
