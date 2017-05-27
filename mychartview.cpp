#include "mychartview.hpp"
#include <QtGui/QMouseEvent>

// this file was taken from Qt5 examples

myChartView::myChartView(QChart *chart, QWidget *parent) :
    QChartView(chart, parent),
    isTouching(false)
{
    setRubberBand(QChartView::RectangleRubberBand);
}

bool myChartView::viewportEvent(QEvent *event)
{
    if (event->type() == QEvent::TouchBegin) {
        isTouching = true;
        chart()->setAnimationOptions(QChart::NoAnimation);
    }
    return QChartView::viewportEvent(event);
}

void myChartView::mousePressEvent(QMouseEvent *event)
{
    if (isTouching)
        return;
    QChartView::mousePressEvent(event);
}

void myChartView::mouseMoveEvent(QMouseEvent *event)
{
    if (isTouching)
        return;
    QChartView::mouseMoveEvent(event);
}

void myChartView::mouseReleaseEvent(QMouseEvent *event)
{
    if (isTouching)
        isTouching = false;

    chart()->setAnimationOptions(QChart::SeriesAnimations);

    QChartView::mouseReleaseEvent(event);
}


void myChartView::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_Plus:
        chart()->zoomIn();
        break;
    case Qt::Key_Minus:
        chart()->zoomOut();
        break;

    case Qt::Key_Left:
        chart()->scroll(-10, 0);
        break;
    case Qt::Key_Right:
        chart()->scroll(10, 0);
        break;
    case Qt::Key_Up:
        chart()->scroll(0, 10);
        break;
    case Qt::Key_Down:
        chart()->scroll(0, -10);
        break;
    default:
        QGraphicsView::keyPressEvent(event);
        break;
    }
}
