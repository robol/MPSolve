import QtQuick 2.1
import QtQuick.Controls 1.0

Row {

    signal solveRequested(string polyText);

    spacing: 16
    width: parent.width
    height: 32

    TextField {
        id: polyInput
        placeholderText: "Insert a polynomial here"
        height: solveButton.height
        width: parent.width - 100
    }

    Button {
        id: solveButton
        text: "Solve"
        width: 84

        onClicked: {
            solveRequested(polyInput.text);
        }
    }

}
