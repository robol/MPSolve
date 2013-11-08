import QtQuick 2.0
import Ubuntu.Components 0.1

Row {

    spacing: units.gu(1)

    TextField {
        id: polyInput
        placeholderText: "Insert polynomial here"
        height: solveButton.height
    }

    Button {
        id: solveButton
        text: "Solve"

        onClicked: {
            solver.solvePoly(polyInput.text);
        }
    }

}
