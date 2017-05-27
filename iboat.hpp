#ifndef IBOAT_HPP
#define IBOAT_HPP

#include <QObject>
#include <QList>
#include <QDateTime>


#include <iboattrajectory.hpp>
#include <iboattrajectorydata.hpp>


class iBoat
{
    public:
        iBoat();


        // object for holding trajecotry data
        IBoatTrajectoryData trajectoryData;

        // only call function for loading trajectory file from source
        bool loadIBoatHurricane(QString file);
        bool loadIBoatPortoTaxi(QString file);

        // start iBoat
        // call function fot inicialization and then start method itself
        void runIBoat();

        // set parameters
        void setParam(double Theta, double Lambda, qint64 CellNumberX, qint64 CellNumberY);

    private:

        // iBoat parameters
        double Theta;
        double Lambda;

        // start detection over all data
        void process();
        // calculate support for actual point
        double support(QList<TrajectoryIBoat *> *hasPathSet, QList<TrajectoryIBoat *> *workingSet);
        // retrun position of first find cell from listNeighborsCells in traj
        qint64 pos(TrajectoryIBoat *traj, QList<QPoint> listNeighborsCells);
        // check if actual point is anomalous
        bool isThetaAnomalous(double supportV);
        // compute sigma for anomalous score
        double sigma(double supportV);
        // return list of neighbor cells of cell including cell
        QList<QPoint> neighborsCells(CellIBoat *cell);
        // return all trajectories in workingSet wich has same path as cells in window
        QList<TrajectoryIBoat *> hasPathWindow(QList<TrajectoryIBoat *> workingSet, QList<CellIBoat *> window);
};

//compute euslidean distance
double euclideanDistanceI(QPointF startPoint, QPointF endPoint);

#endif // IBOAT_HPP
