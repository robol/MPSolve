#ifndef XMPSOLVE_QROOTSRENDERER_H
#define XMPSOLVE_QROOTSRENDERER_H

#include <QWidget>
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

    /**
     * @brief zoomIn handles the zooming operations triggering update().
     */
    void zoomIn();

    /**
     * @brief zoomOut handles the zoomin triggering update().
     */
    void zoomOut();

    /**
     * @brief setCenter handled the recentering calling update()
     * @param x The x coordinate of the new center point.
     * @param y The y coordinate of the new center point.
     */
    void setCenter(double x, double y);

    void mousePressEvent(QMouseEvent * event);
    void mouseReleaseEvent(QMouseEvent * event);
    void mouseMoveEvent(QMouseEvent * event);

private slots:
    void reloadRootsWrapper();

private:
    bool mDragging;
    QPointF mPreviousPosition;

signals:

public slots:

};

} // namespace xmpsolve

#endif // XMPSOLVE_QROOTSRENDERER_H
