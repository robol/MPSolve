#include "rootsrenderer.h"
#include "rootsmodel.h"
#include <QPainter>
#include <QDebug>
#include <cmath>

namespace xmpsolve {

RootsRenderer::RootsRenderer()
{
    m_maxImagModule = m_maxRealModule = 0.0;
    m_model = NULL;
    mCenter = QPointF(0.0, 0.0);
}

void
RootsRenderer::reloadRoots()
{
    m_roots.clear();
    m_maxRealModule = m_maxImagModule = DBL_MIN;

    for (int i = 0; i < m_model->rowCount(); i++)
    {
        Root * root = (Root*) (m_model->data(m_model->index(i), RootsModel::ROOT).value<void*>());

        // Keep track of the maximum module that we reach in both directions, to
        // create a sensible plot later.
        m_maxRealModule = qMax (fabs(root->get_real_part()), m_maxRealModule);
        m_maxImagModule = qMax (fabs(root->get_imag_part()), m_maxImagModule);

        m_roots.append(QPointF(root->get_real_part(), root->get_imag_part()));
    }
}

void
RootsRenderer::zoomIn()
{
    m_maxRealModule = m_maxRealModule / sqrt(2);
    m_maxImagModule = m_maxImagModule / sqrt(2);
}

void
RootsRenderer::zoomOut()
{
    m_maxRealModule = m_maxRealModule * sqrt(2);
    m_maxImagModule = m_maxImagModule * sqrt(2);
}

void
RootsRenderer::setCenter(double x, double y)
{
    mCenter = QPointF(x, y);
}

QPointF
RootsRenderer::scalePoint(QPointF point, int width, int height)
{
    double margin = 24;
    double maxModule = qMax (m_maxRealModule, m_maxImagModule);

    double x = (.5 + .5 * (point.x() - mCenter.x()) / maxModule) * (width - 2 * margin);
    double y = (.5 - .5 * (point.y() - mCenter.y()) / maxModule) * (height - 2 * margin);

    QPointF newPoint(x + margin, y + margin);
    return newPoint;
}

QPointF
RootsRenderer::scaleVectorInverse(QPointF point, int width, int height)
{
    double margin = 24;
    double maxModule = qMax (m_maxRealModule, m_maxImagModule);

    double x = point.x() / (width  - 2 * margin) * maxModule * 2.0;
    double y = - point.y() / (height - 2 * margin) * maxModule * 2.0;

    return QPointF(x, y);
}

QPointF
RootsRenderer::scaleVector(QPointF point, int width, int height)
{
    double margin = 24;
    double maxModule = qMax (m_maxRealModule, m_maxImagModule);

    double x = .5 * point.x() / maxModule * (width - 2 * margin);
    double y = -.5 * point.y() / maxModule * (height - 2 * margin);

    return QPointF(x, y);
}

void
RootsRenderer::drawTicks(QPainter& painter, double w, double h)
{
    double maxModule = qMax(m_maxImagModule, m_maxRealModule);

    if (maxModule == 0.0 || m_roots.length() == 0)
        return;

    int tick_distance_eps = log10 (maxModule * 2);
    double tick_distance = pow(10, tick_distance_eps);

    // Adjust the tick_distance
    if (tick_distance > maxModule / 3)
        tick_distance /= 3.0;

    // Set a small enough font so numbers don't overlap.
    // painter.setFont(QFont(font().family(), 9));
    painter.setPen(QColor(Qt::gray));

    for(int i = -10; i <= 10; i++) {
        // X axis
        QPointF tick_center = scalePoint(QPointF(tick_distance * i, 0), w, h);
        painter.drawLine(tick_center - QPointF(0,4), tick_center + QPointF(0,4));

        QString text;
        text.sprintf("%1.2e", tick_distance * i);
        painter.drawText(QRectF(tick_center - scaleVector(QPointF(tick_distance / 2, 0), w, h) + QPointF(0,-12),
                                tick_center + scaleVector(QPointF(tick_distance / 2, 0), w, h)),
                                Qt::AlignCenter, text);

        // No double zeros, please
        if (i == 0)
            continue;

        // Y axis
        tick_center = scalePoint(QPointF(0, tick_distance * i), w, h);
        painter.drawLine(tick_center - QPointF(4,0), tick_center + QPointF(4,0));

        text.sprintf("%1.2e", tick_distance * i);
        QRectF br = painter.boundingRect(QRectF(tick_center + scaleVector(QPointF(0, tick_distance / 2), w, h),
                                               tick_center + scaleVector(QPointF(0, tick_distance / 2), w, h) + QPointF(50, 0)),
                                        Qt::AlignLeft, text);
        painter.drawText(QRectF(tick_center + scaleVector(QPointF(0, tick_distance / 2), w, h),
                                tick_center + scaleVector(QPointF(0, tick_distance / 2), w, h) +
                                QPointF(br.right() + 4, 0)),
                                Qt::AlignCenter, text);

    }
}

void
RootsRenderer::handlePaintEvent(QPainter& painter, int w, int h, QPaintEvent *)
{
    double axis_margin = 6;

    // Draw the background
    painter.setBrush(QColor::fromRgb(255, 255, 255));
    painter.drawRect(QRect(0, 0, w, h));

    // Draw the cartesian axis
    painter.setPen(QColor(Qt::gray));

    QPointF xMax = scalePoint(QPointF(m_maxImagModule, 0), w, h) - QPointF(axis_margin, 0.0);
    QPointF xMin = scalePoint(QPointF(-m_maxImagModule, 0), w, h) + QPointF(axis_margin, 0.0);
    QPointF yMax = scalePoint(QPointF(0, m_maxImagModule), w, h) - QPointF(0.0, axis_margin);
    QPointF yMin = scalePoint(QPointF(0, -m_maxImagModule), w, h) + QPointF(0.0, axis_margin);

    painter.drawLine(xMin, xMax);
    painter.drawLine(yMin, yMax);

    // Draw the arrows on top of them
    painter.drawLine(xMax - QPointF(8.0, 4.0), xMax);
    painter.drawLine(xMax, xMax - QPointF(8.0, -4.0));
    painter.drawLine(xMin + QPointF(8.0, 4.0), xMin);
    painter.drawLine(xMin, xMin + QPointF(8.0, -4.0));

    painter.drawLine(yMax + QPointF(4.0, 8.0), yMax);
    painter.drawLine(yMax + QPointF(-4.0, 8.0), yMax);
    painter.drawLine(yMin - QPointF(4.0, 8.0), yMin);
    painter.drawLine(yMin - QPointF(-4.0, 8.0), yMin);

    drawTicks(painter, w, h);

    // The default color for the point is red.
    painter.setBrush(QBrush("red"));
    painter.setPen(QColor(Qt::red));

    // ...and we need this to keep track of the MARKED point,
    // if it exists.
    int markedPoint = -1;

    for (int i = 0; i < m_roots.length(); i++)
    {
        QPointF p = scalePoint(m_roots[i], w, h);

        // We draw only not MARKED points. The MARKED one will be drawn in
        // green outside of this for loop.
        if (! m_model->data(m_model->index(i), RootsModel::MARKED).toBool())
            painter.drawEllipse(p, 2, 2);
        else
            markedPoint = i;
    }

    // If we have left out the MARKED point (and this happend only if this
    // point exists) then draw it with green, and bigger than the others.
    if (markedPoint != -1)
    {
        QPointF p = scalePoint(m_roots[markedPoint], w, h);
        painter.setBrush(QBrush("green"));
        painter.setBrush(QColor(Qt::green));
        painter.setPen(QColor(Qt::green));
        painter.drawEllipse(p, 3, 3);
    }
}

} // namespace xmpsolve
