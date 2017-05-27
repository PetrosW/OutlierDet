#ifndef MYCHART_HPP
#define MYCHART_HPP

#include <QtCharts/QChart>

QT_BEGIN_NAMESPACE
class QGestureEvent;
QT_END_NAMESPACE

QT_CHARTS_USE_NAMESPACE

// this file was taken from Qt5 examples

class myChart : public QChart

{
public:
    explicit myChart(QGraphicsItem *parent = 0, Qt::WindowFlags wFlags = 0);
    ~myChart();

protected:
    bool sceneEvent(QEvent *event);

private:
    bool gestureEvent(QGestureEvent *event);

private:

};

#endif // MYCHART_HPP
