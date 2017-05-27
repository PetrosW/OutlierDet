#ifndef IBOATTRAJECTORY_HPP
#define IBOATTRAJECTORY_HPP

#include <QObject>
#include <QPoint>
#include <QList>

class TrajectoryIBoat;

class CellIBoat
{
    public:
        CellIBoat(QPoint cell);
        // X and Y of cell
        QPoint cell;
        // trajectories wich are in cell
        // where key is position of cell in trajectory
        QMap<qint64, QList<TrajectoryIBoat *>> listTrajectory;
        // add new trajecoty to cell
        void addTraj(qint64 pos, TrajectoryIBoat *traj);
};

class TrajectoryIBoat
{
    public:
        TrajectoryIBoat(qint64 ID);
        qint64 ID;
        // list of cell via trajectories go
        QList<CellIBoat *> listCells;
        // list of trajectory points
        QList<QPointF> listPoints;
        // anomali score for every point
        QList<double> score;
        // list of anomalous points
        QList<QPointF> listAnomalousPoints;
        // list of index to listPoints
        // for anomalous points
        QList<qint64> listAnomalousIndex;
        // time of trajectory start
        double time;
};



#endif // IBOATTRAJECTORY_HPP
