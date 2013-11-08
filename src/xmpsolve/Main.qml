import QtQuick 2.0
import Ubuntu.Components 0.1
import Ubuntu.Components.ListItems 0.1

MainView {
    id: root

    objectName: "mainView"
    applicationName: "MPSolve"

    width: units.gu(50)
    height: units.gu(75)

    Page {
        id: polyInputPage
        title: "MPSolve"

        Column {
            id: pageLayout

            anchors.fill: parent
            anchors.margins: units.gu(1)
            spacing: units.gu(2)

            PolyInputField {}

            Text { text: "Approximations computed by MPSolve: " }

            ApproximationList {
                model: rootsModel
            }


        }

    }
}
