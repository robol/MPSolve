import QtQuick 2.0

Canvas {

    id: canvas
    height: 400
    width: 400

    function drawBorder(ctx) {
        ctx.moveTo(0,0)
        ctx.lineTo(canvas.width, 0);
        ctx.lineTo(canvas.width, canvas.height);
        ctx.lineTo(0, canvas.height);
        ctx.lineTo(0, 0);
        ctx.stroke();
    }

    function coordsToPoint(x, y) {

    }

    property int prevX
    property int prevY
    property int lineWidth: 2
    property color drawColor: "black"

    onPaint: {
        var ctx = getContext('2d');
        drawBorder(ctx);

        var nroots = rootsModel.rowCount();
        var i = 0;
        for (i = 0; i < nroots; i++) {
            var x = rootsModel.getPointX(i);
            var y = rootsModel.getPointY(i);


        }
    }

    onCanvasSizeChanged: {
        requestPaint();
    }

}
