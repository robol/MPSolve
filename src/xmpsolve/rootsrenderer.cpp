#include "rootsrenderer.h"
#include <QPainter>
#include <QDebug>
#include <cmath>

namespace xmpsolve {

RootsRenderer::RootsRenderer(QWidget *parent) :
    QWidget(parent)
{
    m_maxImagModule = m_maxRealModule = 0.0;
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
    double maxModule = qMax (m_maxRealModule, m_maxImagModule);
    QPointF center(width / 2, height / 2);

    double x = (.5 + .5 * point.x() / maxModule) * (width - 2 * margin);
    double y = (.5 - .5 * point.y() / maxModule) * (height - 2 * margin);

    QPointF newPoint(x + margin, y + margin);
    return newPoint;
}

void
RootsRenderer::drawTicks(QPainter& painter)
{
    double w = width();
    double h = height();

    double maxModule = qMax(m_maxImagModule, m_maxRealModule);

    if (maxModule == 0.0)
        return;

    int tick_distance_eps = log10 (maxModule);
    double tick_distance = pow(10, tick_distance_eps);

    // Adjust the tick_distance
    if (tick_distance > maxModule / 3)
        tick_distance /= 3.0;

    // Set a small enough font so numbers don't overlap.
    painter.setFont(QFont(font().family(), 9));
    painter.setPen(QColor(Qt::gray));

    for(int i = -10; i <= 10; i++) {
        // X axis
        QPointF tick_center = scalePoint(QPointF(tick_distance * i, 0), w, h);
        painter.drawLine(tick_center - QPointF(0,4), tick_center + QPointF(0,4));

        QString text;
        text.sprintf("%1.2e", tick_distance * i);
        painter.drawText(QRectF(tick_center - scalePoint(QPointF(tick_distance / 2, 0), w, h) + QPointF(0,-12),
                                tick_center + scalePoint(QPointF(tick_distance / 2, 0), w, h)),
                                Qt::AlignCenter, text);

        // No double zeros, please
        if (i == 0)
            continue;

        // Y axis
        tick_center = scalePoint(QPointF(0, tick_distance * i), w, h);
        painter.drawLine(tick_center - QPointF(4,0), tick_center + QPointF(4,0));

        text.sprintf("%1.2e", tick_distance * i);
        QRectF br = painter.boundingRect(QRectF(tick_center - scalePoint(QPointF(0, tick_distance / 2), w, h),
                                               tick_center + scalePoint(QPointF(0, tick_distance / 2), w, h) + QPointF(50, 0)),
                                        Qt::AlignLeft, text);
        painter.drawText(QRectF(tick_center - scalePoint(QPointF(0, tick_distance / 2), w, h),
                                tick_center + scalePoint(QPointF(0, tick_distance / 2), w, h) +
                                QPointF(br.right() + 4, 0)),
                                Qt::AlignCenter, text);

    }
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
    painter.setPen(QColor(Qt::gray));
    painter.drawLine(QPointF(axis_margin, h/2), QPointF(w - axis_margin, h/2));
    painter.drawLine(QPointF(w/2, axis_margin), QPointF(w/2, h - axis_margin));

    // Draw the arrows on top of them
    painter.drawLine(QPointF(w - axis_margin, h/2), QPointF(w - 8 - axis_margin, h/2 - 4));
    painter.drawLine(QPointF(w - axis_margin, h/2), QPointF(w - 8 - axis_margin, h/2 + 4));
    painter.drawLine(QPointF(w/2, axis_margin), QPointF(w/2 + 4, 8 + axis_margin));
    painter.drawLine(QPointF(w/2, axis_margin), QPointF(w/2 - 4, 8 + axis_margin));

    drawTicks(painter);

    painter.setBrush(QBrush("red"));
    painter.setPen(QColor(Qt::red));
    foreach(QPointF point, m_roots)
    {
        QPointF p = scalePoint(point, w, h);
        painter.drawEllipse(p, 2, 2);
    }
}

} // namespace xmpsolve
