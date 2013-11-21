import QtQuick 2.1
import QtQuick.Controls 1.0

ApplicationWindow {
    id: root

    title: "MPSolve"

    width: 400
    height: 600

    visible: true

    Column {
        id: pageLayout

        anchors.fill: parent
        anchors.margins: 16
        spacing: 32

        Text  {
            text: PACKAGE_STRING
            font.pointSize: 18
        }

        PolyInputField {}

        ApproximationList {
            model: rootsModel
        }


    }
}
