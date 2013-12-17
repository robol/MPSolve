import QtQuick 2.1
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0

ApplicationWindow {
    id: root
    title: PACKAGE_STRING

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
            font.pointSize: 16
        }

        PolyInputField {}

        // Phone view: activated only if root.width < root.height
        SwitchableApproximationView {
            id: switchableView
            width: parent.width
            height: root.height - 175
            visible: (root.width < root.height)
        }

        // Tablet view: activated only if root.width >= root.height
        CombinedApproximationView {
            width: root.width
            height: root.height - 140
            visible: (root.width >= root.height)
        }
    }




}
