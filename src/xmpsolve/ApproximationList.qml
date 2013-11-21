import QtQuick 2.1

Column {

    /**
     * @brief Update the description of the approximations.
     *
     * @param count The number of items in the list.
     */
    function updateText(count) {
       if (count > 0)
           approximationTitle.text = "List of approximations: ";
       else
           approximationTitle.text = "No approximations computed."
    }

    property alias model : resultsView.model
    spacing: 16

    Text {
        id: approximationTitle
        text: "No approximations computed."
    }

    ListView {
        id: resultsView
        delegate: approximationDelegate
        spacing: 12

        height: 720
        width: parent.width

        onCountChanged: {
            updateText(count)
        }

        Component {
            id: approximationDelegate

            Item {
                width: parent.width
                height: approximationItem.height

                Text {
                    id: approximationItem
                    text: short_approximation
                    height: 8
                }
            }
        }
    }
}
