import QtQuick 2.0
import QtQuick.Controls 1.0
import MPSolve 1.0

Column {

    function toggleView() {
        switch (stackView.state) {
            case 0:
                stackView.push(approximationList)
                break;
            case 1:
                stackView.push(rootsRenderer);
                break;
        }

        stackView.state = (stackView.state + 1) % 2;
    }

    spacing: 16

    StackView {
        id: stackView
        initialItem: rootsRenderer

        property int state: 0;

        ApproximationList {
            id: approximationList
            model: rootsModel
            clip: true
            visible: false
        }

        QQuickRootsRenderer {
            id: rootsRenderer
            model: rootsModel
        }
    }

    Button {
        text: "Toggle Graphics / Approximation list"
        width: parent.width
        height: 32

        onClicked: toggleView()
    }

}
