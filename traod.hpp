#ifndef TRAOD_HPP
#define TRAOD_HPP

#include <QObject>
#include <QPoint>
#include <QtMath>
#include <QVector>
#include <QVector2D>
#include <QDateTime>

#include <QDebug>

class CoarseSegment;
class Segment;
class TRAODTrajecotryData;
#include <traodtrajectory.hpp>
#include <traodtrajecotrydata.hpp>

class TRAOD
{

    public:

        TRAOD();
        // start TRAOD
        // call function fot inicialization and then start method itself
        void runTRAOD();
        // only call function for loading trajectory file from source
        bool loadTRAODHurricane(QString file);
        bool loadTRAODPortoTaxi(QString file);
        // set parameters
        void setParams(double p, double F, double D, double wa, double wper, double wpar);
        TRAODTrajecotryData trajectoryData;


    private:
        double MDLCostAdjustment;
        double p;
        double F;
        double D;

        double wa;
        double wper;
        double wpar;

        // standard deviation of distance between fine segments
        double stdDev;
        // average count of fine segment within of standard deviation distance
        double averageDensity;
        double sumDensity;
        qint64 numOfAllSegments;
        qint64 numOfAllCoarseSegments;

        // global lists for close fine segment of every trajectory against other one trajectory
        QList<Segment *> *cl1;
        QList<QList<Segment *> *> *cl2;

        // start detection over all data with optimization version of TRAOD
        void processOrigOpt();
        // find all fine segment of trajIn of coarse segment within their bounds
        void processSegmentOrigOpt(CoarseSegment *coarseSegment, TrajectoryTRAOD *trajIn);
        // return lower and upper bound between coarseSeg1 and coarseSeg2
        QPair<double, double> getBounds(CoarseSegment *coarseSeg1, CoarseSegment *coarseSeg2);
        // return list of close fine segment for every trajectory
        QList<qint64> getNumOfCloseTOpt();
        double computeAdjustingCoefficientOpt(Segment *seg);
        bool isTrajectoryOutlier(QList<Segment *> *listSegments);

        // return projection point in (LS, LE) segment from LP point
        QPointF perpenPoint(QPointF LS, QPointF LE, QPointF LP, double* k);
        // return final distance
        double disAllFin(Segment *L1, Segment *L2);
        // return angle distance as defined in TRAOD, through angle parameter return angel between segment
        double disAngle(double L2, QPointF L1S, QPointF L1E, QPointF L2S, QPointF L2E, double *angle);
        // return perpendicular distance as defined in TRAOD
        double disPerpen(double lk1, double lk2);
        // return paraller distance from p1 and p2 as defined in TRAOD
        QPair<double, double> disPar(double k1, double k2, QPointF S, QPointF E, QPointF p1, QPointF p2);
        // functions for MDL costs
        double MDLnopar(Segment *L1, Segment *L2);
        // coarseSegParam return new coarse segment properties, later used for bounds
        double MDLpar(qint64 startIndex, qint64 endIndex, QList<Segment *> *listSegments, QVector<double> &coarseSegParam);
        // apply weigth parameters to distances
        double weightDis(double disAngle, double disPer, double disPar);
        double standardDeviation();
        // approximation function for standard deviation
        double standardDeviationAprox();
        void densityForAllSegmenst();
        // add coarse segments to trajectories from fine segments
        void createCoarseSegments();
};

double euclideanDistance(QPointF startPoint, QPointF endPoint);

#endif // TRAOD_HPP
