import QtQuick 2.0
import QtQuick.Layouts 1.0

ColumnLayout {

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
