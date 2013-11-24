#include "mainqmlview.h"
#include "qquickrootsrenderer.h"
#include <QQmlEngine>
#include <QtQml>

namespace xmpsolve {

MainQmlView::MainQmlView()
{
    m_rootContext = rootContext();
    m_model = m_solver.rootsModel();
    inflateObjects();
    loadTypes();

    load(QUrl("qrc:/qml/Main.qml"));
}

void
MainQmlView::loadTypes()
{
    qmlRegisterType<QQuickRootsRenderer>("MPSolve", 1, 0, "QQuickRootsRenderer");
}

void
MainQmlView::inflateObjects()
{
    QQmlContext *ctx = m_rootContext;

    ctx->setContextProperty("rootsModel", m_model);
    ctx->setContextProperty("solver", &m_solver);
    ctx->setContextProperty("parser", &m_parser);

    // Set some other cosmetic values
    ctx->setContextProperty("PACKAGE_STRING", PACKAGE_STRING);
}

} // namespace xmpsolve
