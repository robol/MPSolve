#ifndef XMPSOLVE_MAINQMLVIEW_H
#define XMPSOLVE_MAINQMLVIEW_H

#include <QtQml/QQmlApplicationEngine>
#include <QQmlContext>
#include "rootsmodel.h"
#include "rootsrenderer.h"
#include "polynomialsolver.h"
#include "polynomialparser.h"

namespace xmpsolve {

class MainQmlView : public QQmlApplicationEngine
{
    Q_OBJECT
public:
    explicit MainQmlView();

private:
    /**
     * @brief m_model is the model that holds the current approximations
     * computed by MPSolve.
     */
    RootsModel *m_model;

    /**
     * @brief m_solver is the default polynomial solver that can be used
     * to perform the computations.
     */
    PolynomialSolver m_solver;

    /**
     * @brief m_parser is a parse for polynomials that transforms strings into Polynomial
     * instances. Those can then be fed into m_solver.
     */
    PolynomialParser m_parser;

    /**
     * @brief m_rootContext is the root Context of the current QQmlEngine
     * embedded in the QQuickView.
     *
     * It is used in here to allow to the Qml Items an easy access to
     * C++ objects.
     */
    QQmlContext *m_rootContext;

    /**
     * @brief inflateObjects loads the useful members of the View in the root Context
     * of the QQmlEngine.
     */
    void inflateObjects();

signals:

public slots:

};

} // namespace xmpsolve

#endif // XMPSOLVE_MAINQMLVIEW_H
