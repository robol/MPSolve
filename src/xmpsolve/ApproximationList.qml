import QtQuick 2.0

ListView {
    id: resultsView
    delegate: approximationDelegate
    spacing: units.gu(2)

    height: units.gu(45)
    width: parent.width

    Component {
        id: approximationDelegate

        Item {
            width: parent.width
            height: approximationItem.height

            Text {
                id: approximationItem
                text: short_approximation

                height: units.gu(2)
            }
        }
    }
}
