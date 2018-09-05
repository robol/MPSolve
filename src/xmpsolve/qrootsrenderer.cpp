#include "qrootsrenderer.h"
#include <QDebug>
#include <QMouseEvent>

namespace xmpsolve {

QRootsRenderer::QRootsRenderer(QWidget *parent) :
    QWidget(parent)
{
    mDragging = false;
}

void
QRootsRenderer::setModel(RootsModel *model)
{
    m_model = model;
    reloadRoots();

    // Grab data changes and model resets.
    m_model->connect(m_model, SIGNAL(dataChanged(QModelIndex,QModelIndex)),
                     this, SLOT(reloadRootsWrapper()));
    m_model->connect(m_model, SIGNAL(modelReset()), this, SLOT(reloadRootsWrapper()));
}

void
QRootsRenderer::reloadRootsWrapper()
{
    reloadRoots();
    update();
}

void
QRootsRenderer::zoomIn()
{
    RootsRenderer::zoomIn();
    update();
}

void
QRootsRenderer::zoomOut()
{
    RootsRenderer::zoomOut();
    update();
}

void
QRootsRenderer::setCenter(double x, double y)
{
    RootsRenderer::setCenter(x, y);
    update();
}

void
QRootsRenderer::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        mDragging = true;
        mPreviousPosition = event->localPos();
    }
    QWidget::mousePressEvent(event);
}

void
QRootsRenderer::mouseReleaseEvent(QMouseEvent *event)
{
    if (mDragging) {
        mDragging = false;
    }
    QWidget::mouseReleaseEvent(event);
}

void
QRootsRenderer::mouseMoveEvent(QMouseEvent *event)
{
    if (mDragging) {
        QPointF distance = event->localPos() - mPreviousPosition;
        distance = scaleVectorInverse(distance, width(), height());
        setCenter(center().x() - distance.x(), center().y() - distance.y());
        mPreviousPosition = event->localPos();
    }
    QWidget::mouseMoveEvent(event);
}

void
QRootsRenderer::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);    
    painter.setRenderHint(QPainter::Antialiasing, true);
    handlePaintEvent(painter, width(), height(), event);
}

} // namespace xmpsolve
