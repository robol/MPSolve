import QtQuick 2.1
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0
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

            width: parent.width
            height: parent.height - 150

            id: viewApproximationRow
            spacing: 12

            ApproximationList {
                id: approximationList
                model: rootsModel
                clip: true

                width: 200
                height: parent.height
            }

            QQuickRootsRenderer {
                model: rootsModel

                width: parent.width - 212
                height: parent.height
            }

        }


    }
}
