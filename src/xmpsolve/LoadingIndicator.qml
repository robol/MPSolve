import QtQuick 2.1
import QtQuick.Layouts 1.0

Rectangle {

    color: "#aaeaaa"

    width: parent.width
    height: parent.height

    ColumnLayout {

        anchors.fill: parent

        Image {
            id: loading
            source: "../images/loading.png"
            fillMode: Image.PreserveAspectFit

            Layout.alignment: Qt.AlignHCenter | Qt.AlignVCenter

            onVisibleChanged: {
                loading.width = 128;
                loading.height = 128;
            }

            NumberAnimation on rotation {
                 from: 0
                 to: 360
                 running: loading.visible == true
                 loops: Animation.Infinite
                 duration: 900
             }
        }

        Text {
            text: "MPSolve is working for you ..."
            font.pointSize: 16
            Layout.alignment: Qt.AlignHCenter
        }

    }
}
