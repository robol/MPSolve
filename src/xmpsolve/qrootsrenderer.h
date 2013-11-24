#ifndef XMPSOLVE_QROOTSRENDERER_H
#define XMPSOLVE_QROOTSRENDERER_H

#include <QWidget>
#include <QPainter>
#include "rootsrenderer.h"

namespace xmpsolve {

class QRootsRenderer : public QWidget, public RootsRenderer
{
    Q_OBJECT
public:
    explicit QRootsRenderer(QWidget *parent = 0);

    /**
     * @brief setRoots can be used to set the roots that the RootsRenderer
     * shall renderer.
     * @param model The roots model that should be displayed.
     */
    void setModel(RootsModel* model);

    void paintEvent(QPaintEvent* event);

private slots:
    void reloadRootsWrapper();

signals:

public slots:

};

} // namespace xmpsolve

#endif // XMPSOLVE_QROOTSRENDERER_H
