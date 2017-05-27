#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>
#include <QMessageBox>

#include <QDebug>

#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QChart>
#include <QtCharts/QScatterSeries>
#include <QTableWidget>
#include <QTableWidgetItem>

#include <mychart.hpp>
#include <mychartview.hpp>
#include <traod.hpp>
#include <iboat.hpp>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);

    // objects of detection methods
    TRAOD traodObject;
    iBoat iBoatObject;

    // list of series for graphs for selected anomalous trajectory
    // in iBoat method
    QList<QtCharts::QScatterSeries *> listAnomalySeries;

    // method for printing graphs for TRAOD method
    void printTrajectoriesTRAOD();
    void printCoarseTrajectoryTRAOD();
    void printTrajectoryOutlierTRAOD();
    void printSegmentOutlierTRAOD();

    // object of graphs
    myChart *chartTRAOD;
    myChart *chartIBoat;
    myChart *chartIBoatScore;

    // method for printing graphs for iBoat method
    void printTrajIBoat();
    void printAnomalousTrajIBoat();
    void printCellIBoat();
    void printAnomalousPointIBoat();
    void fillTableList();

    void setParamIBoat();
    void setParamTRAOD();

    ~MainWindow();

public slots:
    // slots for every method with different data set
    // if you need new data set copy one of these method
    // and change data load function to your own new
    void runTRAODHurricane();
    void runTRAODPortoTaxi();
    void runIBoatHurricane();
    void runIBoatPortoTaxi();

    // when select anomalous trajectory from table (iBoat)
    // show it in graph and hide previous trajectory
    void cellIBoatSelected(qint64 row, qint64 column);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_HPP
