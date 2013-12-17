import QtQuick 2.0
import MPSolve 1.0

Row {

    spacing: 16

    ApproximationList {
        model: rootsModel
        width: 350
        height: parent.height
    }

    QQuickRootsRenderer {
        model: rootsModel
        width: root.width - 400
        height: parent.height
    }
}
