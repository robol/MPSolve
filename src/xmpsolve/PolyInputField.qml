import QtQuick 2.1
import QtQuick.Controls 1.0

Column {

    spacing: 8

    Row {

        spacing: 16

        TextField {
            id: polyInput
            placeholderText: "Insert a polynomial here"
            height: solveButton.height
            width: 500
        }

        Button {
            id: solveButton
            text: "Solve"

            onClicked: {
                solver.solvePoly(polyInput.text);
            }
        }

    }

}
