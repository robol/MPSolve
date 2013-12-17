import QtQuick 2.1
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0

ApplicationWindow {
    id: root
    title: PACKAGE_STRING

    width: 480
    height: 600

    visible: true
    color: "#cacafa"

    // Hook up connections to the solver, so we can handle the spinning
    // bar and other small thinkgs.
    Connections {
        target: solver

        onSolved: {
            loadingIndicator.visible = false;
            mainStack.push(mainView);
        }
    }

    StackView {

        id: mainStack

        MainView {
            id: mainView

            onSolveRequested: {
                mainStack.push(loadingIndicator);
                mainView.visible = false;
                solver.solvePoly(polyText);
            }
        }

        LoadingIndicator {
            id: loadingIndicator
            visible: false
        }
    }


}
