#include "mainqmlview.h"
#include <QQmlEngine>
#include <QtQml>

namespace xmpsolve {

MainQmlView::MainQmlView()
{
    m_rootContext = engine()->rootContext();
    m_model = m_solver.rootsModel();
    inflateObjects();

    // Load the Qml file after having inflated the objects, otherwhile
    // we will get a lot of undefined references errors.
    setSource(QUrl("qrc:/qml/Main.qml"));
}

void
MainQmlView::inflateObjects()
{
    QQmlContext *ctx = m_rootContext;

    ctx->setContextProperty("rootsModel", m_model);
    ctx->setContextProperty("solver", &m_solver);
    ctx->setContextProperty("parser", &m_parser);
}

} // namespace xmpsolve
