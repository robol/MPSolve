#include "rootsrenderer.h"
#include <QPainter>
#include <QDebug>
#include <cmath>

namespace xmpsolve {

RootsRenderer::RootsRenderer(QWidget *parent) :
    QWidget(parent)
{
}

void
RootsRenderer::setRoots(QList<Root *> roots)
{
    m_roots.clear();
    m_maxRealModule = m_maxImagModule = DBL_MIN;

    foreach(Root *root, roots)
    {
        // Keep track of the maximum module that we reach in both directions, to
        // create a sensible plot later.
        m_maxRealModule = qMax (fabs(root->get_real_part()), m_maxRealModule);
        m_maxImagModule = qMax (fabs(root->get_imag_part()), m_maxImagModule);

        m_roots.append(QPointF(root->get_real_part(), root->get_imag_part()));
    }

    update();
}

QPointF
RootsRenderer::scalePoint(QPointF point, int width, int height)
{
    double margin = 24;
    QPointF center(width / 2, height / 2);

    double x = (.5 + .5 * point.x() / m_maxRealModule) * (width - 2 * margin);
    double y = (.5 - .5 * point.y() / m_maxImagModule) * (height - 2 * margin);

    QPointF newPoint(x + margin, y + margin);
    return newPoint;
}

void
RootsRenderer::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    double axis_margin = 6;

    double w = width();
    double h = height();

    // Draw the background
    painter.setBrush(QColor::fromRgb(255, 255, 255));
    painter.drawRect(QRect(0, 0, w, h));

    // Draw the cartesian axis
    painter.setPen(QColor::fromRgb(0,0,0));
    painter.drawLine(QPointF(axis_margin, h/2), QPointF(w - axis_margin, h/2));
    painter.drawLine(QPointF(w/2, axis_margin), QPointF(w/2, h - axis_margin));

    // Draw the arrows on top of them
    painter.drawLine(QPointF(w - axis_margin, h/2), QPointF(w - 8 - axis_margin, h/2 - 4));
    painter.drawLine(QPointF(w - axis_margin, h/2), QPointF(w - 8 - axis_margin, h/2 + 4));
    painter.drawLine(QPointF(w/2, axis_margin), QPointF(w/2 + 4, 8 + axis_margin));
    painter.drawLine(QPointF(w/2, axis_margin), QPointF(w/2 - 4, 8 + axis_margin));

    painter.setBrush(QBrush("red"));
    foreach(QPointF point, m_roots)
    {
        QPointF p = scalePoint(point, w, h);
        painter.drawEllipse(p, 2, 2);
    }
}

} // namespace xmpsolve
