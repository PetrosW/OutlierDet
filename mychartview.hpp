#ifndef MYCHARTVIEW_HPP
#define MYCHARTVIEW_HPP

#include <QtCharts/QChartView>
#include <QtWidgets/QRubberBand>

QT_CHARTS_USE_NAMESPACE

// this file was taken from Qt5 examples

class myChartView : public QChartView
{
public:
    myChartView(QChart *chart, QWidget *parent = 0);

protected:
    bool viewportEvent(QEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);

private:
    bool isTouching;
};

#endif // MYCHARTVIEW_HPP
