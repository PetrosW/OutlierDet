#ifndef TRAODTRAJECTORY_HPP
#define TRAODTRAJECTORY_HPP

#include <QObject>
#include <QPoint>
#include <QDebug>

class Segment
{
    public:
        Segment(QPointF startPoint, QPointF endPoint);

        QPointF startPoint;
        QPointF endPoint;
        double length;
        bool isOutlier;
        qint64 density;
};

class CoarseSegment
{
    public:
        CoarseSegment();

        // list of fine segments in coarse segment
        QList<Segment *> listFineSegments;
        QPointF startPoint;
        QPointF endPoint;
        double length;
        bool isOutlier;
        // coarse segment properties
        // its used for bounds computing
        double maxAngle;
        double maxLength;
        double minLength;
        double maxLk;

        double calculateLength();
};

class TrajectoryTRAOD
{

    public:
        TrajectoryTRAOD();

        double time;
        double score;
        bool isOutlier;

        // points (listPoints), fine segments (listSegments),
        // coarse segments (listCoarseSegments) and coarse points (listCoarsePoints)
        // of trajectory
        QList<QPointF> listPoints;
        QList<Segment *> listSegments;
        QList<QPointF> listCoarsePoints;
        QList<CoarseSegment *> listCoarseSegments;

        void addPoint(QPointF point);
        // from points create fine segments
        qint64 createSegments();
};

double euclideanDistance(QPointF startPoint, QPointF endPoint);


#endif // TRAODTRAJECTORY_HPP
