#include "qrootsrenderer.h"

namespace xmpsolve {

QRootsRenderer::QRootsRenderer(QWidget *parent) :
    QWidget(parent)
{
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
QRootsRenderer::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);    
    painter.setRenderHint(QPainter::Antialiasing, true);
    handlePaintEvent(painter, width(), height(), event);
}

} // namespace xmpsolve
