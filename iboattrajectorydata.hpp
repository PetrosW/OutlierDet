#ifndef IBOATTRAJECTORYDATA_HPP
#define IBOATTRAJECTORYDATA_HPP

#include <QObject>
#include <QtMath>
#include <QMap>
#include <QFile>
#include <QPoint>

#include <QDebug>

#include <iboattrajectory.hpp>

class IBoatTrajectoryData
{
    public:
        IBoatTrajectoryData();

        // structures for holding data for detection
        QList<TrajectoryIBoat *> listTrajectories;
        QMap<QPair<qint64, qint64>, CellIBoat *> listCells;

        // clean lists so detection can start again with different parameters
        void cleanLists();
        // loading data (files) to trajectory structures
        // if you need load your own data create your own function
        // your new function must fill listPoints in TrajectoryIBoat object for every trajectory
        // and object of every trajectory add to listTrajectories
        bool readHurricaneFromFile(QString filePath);
        bool readPortoTaxiFromFile(QString filePath);
        // setter for cell count
        void setSizeOfCell(qint64 CellNumberX, qint64 CellNumberY);
        // fill listCells with objects of cell
        void initCell();

        // check if cells are neighbors
        bool isCellsNeighbor(QPoint cell1, QPoint cell2);
        // compute cell index from trajectories points (GPS or whatever)
        qint64 computeIndexY(double Y);
        qint64 computeIndexX(double X);


    private:
        double maxX, minX, maxY, minY;
        double CellNumberX, CellNumberY;
        double sizeOfCellX, sizeOfCellY;

        // map trajectories to cell
        void mapToCell();
        // compute augumented cell from point to point
        QList<QPoint> augumentCell(QPoint S, QPoint E);

        // compute size of cell from min and max of X and Y
        void computeSizeOfCell();
        // find min and max of X and Y
        void getMinMax();
};



#endif // IBOATTRAJECTORYDATA_HPP
