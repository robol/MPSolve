import QtQuick 2.0

Item {

    id: top

    width: parent.width
    height: parent.height

    signal solveRequested(string polyText);

    Column {
        id: pageLayout

        spacing: 16

        anchors.fill: parent
        anchors.margins: 16

        Text  {
            text: PACKAGE_STRING
            font.pointSize: 16
        }

        PolyInputField {
            onSolveRequested: {
                top.solveRequested(polyText);
            }
        }

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
