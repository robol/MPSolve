import QtQuick 2.1
import QtQuick.Controls 1.0
import MPSolve 1.0

ApplicationWindow {
    id: root

    title: "MPSolve"

    width: 700
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

        Row {

            id: viewApproximationRow
            spacing: 12

            ApproximationList {
                id: approximationList
                model: rootsModel

                width: 200
            }

            QQuickRootsRenderer {
                model: rootsModel

                width: 450
                height: 450
            }

        }


    }
}
