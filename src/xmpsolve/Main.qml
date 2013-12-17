import QtQuick 2.1
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0
import MPSolve 1.0

ApplicationWindow {
    id: root

    title: "MPSolve"

    width: 480
    height: 600

    visible: true

    Column {
        id: pageLayout

        anchors.fill: parent
        anchors.margins: 16
        spacing: 16

        Text  {
            text: PACKAGE_STRING
            font.pointSize: 18
        }

        PolyInputField {}

        SwitchableApproximationView {
            width: parent.width
            height: root.height - 175
        }

    }




}
