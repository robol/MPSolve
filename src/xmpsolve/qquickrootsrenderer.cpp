#include "qquickrootsrenderer.h"

using namespace xmpsolve;

QQuickRootsRenderer::QQuickRootsRenderer(QQuickItem *parent) :
    QQuickPaintedItem(parent)
{
}

void
QQuickRootsRenderer::paint(QPainter * painter)
{
    handlePaintEvent(*painter, width(), height(), NULL);
}

void
QQuickRootsRenderer::setModel(QVariant model)
{
    m_model = model.value<RootsModel*>();

    // Grab data changes and model resets.
    m_model->connect(m_model, SIGNAL(dataChanged(QModelIndex,QModelIndex)),
                     this, SLOT(reloadRootsWrapper()));
    m_model->connect(m_model, SIGNAL(modelReset()), this, SLOT(reloadRootsWrapper()));
}

void
QQuickRootsRenderer::reloadRootsWrapper()
{
    reloadRoots();
    update();
}
