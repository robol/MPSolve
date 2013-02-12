#ifndef xmpsolve_ROOTSRENDERER_H
#define xmpsolve_ROOTSRENDERER_H

#include <QWidget>
#include "root.h"

namespace xmpsolve {

class RootsRenderer : public QWidget
{
    Q_OBJECT
public:
    explicit RootsRenderer(QWidget *parent = 0);

    /**
     * @brief setRoots can be used to set the roots that the RootsRenderer
     * shall renderer.
     * @param roots The roots that should be shown. Can be an  empty list
     * to clear a previous plot.
     */
    void setRoots(QList<Root*> roots);
    
signals:
    
public slots:
    void paintEvent(QPaintEvent *);

private:
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

    QList<QPointF> m_roots;

    double m_maxRealModule;
    double m_maxImagModule;
    
};

} // namespace xmpsolve

#endif // xmpsolve_ROOTSRENDERER_H
