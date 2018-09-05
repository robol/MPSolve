#ifndef xmpsolve_ROOTSRENDERER_H
#define xmpsolve_ROOTSRENDERER_H

#include "root.h"
#include "rootsmodel.h"
#include <QPainter>
#include <QPaintEvent>

namespace xmpsolve {

class RootsRenderer
{

public:
    explicit RootsRenderer();

    void handlePaintEvent(QPainter& painter, int w, int h, QPaintEvent *);

    /**
     * @brief zoomIn changes the zoom to magnify more.
     */
    void zoomIn();

    /**
     * @brief zoomOut changes the zoom to see a wider range.
     */
    void zoomOut();

    /**
     * @brief setCenter selects the current centered point.
     * @param x Horizontal coordinate of the center (aka real part)
     * @param y Vertical coordinate of the center (aka imaginary part).
     */
    void setCenter(double x, double y);

    /**
     * @brief center returns the current center point of the plot.
     * @return A QPointF object representing the real and imaginary part of the
     * current center.
     */
    QPointF center() { return mCenter; }

protected:

    /**
     * @brief reloadRoots reloads the roots from the model.
     */
    void reloadRoots();

    /**
     * @brief scalePoint is used internally to scale, flip and translate a point
     * in a such a way that is plotted properly on the coordinate system on the screen.
     *
     * @param point Is the original point that shall be plotted.
     * @param width Is the current width of the widget.
     * @param height Is the current height of the widget.
     * @return The QPointF that can be plotted.
     */
    QPointF scalePoint(QPointF point, int width, int height);

    QPointF scaleVector(QPointF point, int width, int height);
    QPointF scaleVectorInverse(QPointF point, int width, int height);

    /**
     * @brief drawTicks is used internally to draw ticks on the axis.
     */
    void drawTicks(QPainter& painter, double w, double h);

    /**
     * @brief Points that should be displayed.
     */
    QList<QPointF> m_roots;

    /**
     * @brief m_maxRealModule is the maximum module of the real parts of the roots.
     */
    double m_maxRealModule;

    /**
     * @brief m_maxImagModule is the maximum module of the imaginary parts of the roots.
     */
    double m_maxImagModule;

    /**
     * @brief m_model is the model containing the roots that should be displayed.
     */
    RootsModel* m_model;

private:
    QPointF mCenter;
    
};

} // namespace xmpsolve

#endif // xmpsolve_ROOTSRENDERER_H
