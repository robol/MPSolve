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
        spacing: 12
        clip: true

        height: 420
        width: parent.width

        onCountChanged: {
            updateText(count)
        }

        Component {
            id: approximationDelegate

            Item {
                width: parent.width
                height: approximationItem.height

                Row {
                    spacing: 6

                    Text {
                        id: approximationItem
                        text: short_approximation
                        height: 8
                        width: resultsView.width - 150
                    }

                    Text {
                        text: "Radius: " + radius
                        height: 8
                        width: 64
                    }

                }
            }
        }
    }
}
