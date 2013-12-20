import QtQuick 2.1

Column {

    /**
     * @brief Update the description of the approximations.
     *
     * @param count The number of items in the list.
     */
    function updateText(count) {
       if (count > 0) {
           approximationTitle.text = "List of approximations: ";
           approximationTitle.font.pointSize = "14";
       }
       else {
           approximationTitle.text = "No approximations computed.\nInsert a polynomial to see some data here!"
           approximationTitle.font.pointSize = "10";
       }
    }

    property alias model : resultsView.model
    spacing: 16

    Text {
        id: approximationTitle
    }

    ListView {
        id: resultsView
        delegate: approximationDelegate
        spacing: 6
        clip: true

        height: 420
        width: parent.width

        onCountChanged: {
            updateText(count)
        }

        Component {
            id: approximationDelegate

            Rectangle {

                width: parent.width
                height: 48
                color: "#4d4"
                opacity: .9
                radius: 6

                Column {

                    spacing: 6

                    anchors.margins: 6
                    anchors.fill: parent

                    Text {
                        id: approximationItem
                        text: short_approximation
                        height: 14
                        font.pointSize: 11
                    }

                    Text {
                        text: "Inclusion radius: " + inclusion_radius + ", Status: " + status
                        font.pointSize: 9
                        color: "#444"
                        height: 14
                    }

                }

                MouseArea {
                    anchors.fill: parent
                    onClicked: {
                        rootsModel.markRoot(index);
                        switchableView.toggleView();
                    }
                }
            }
        }
    }
}
