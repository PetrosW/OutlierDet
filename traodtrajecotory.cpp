#include "traodtrajectory.hpp"

Segment::Segment(QPointF startPoint, QPointF endPoint) : startPoint(startPoint), endPoint(endPoint) {
    length = euclideanDistance(startPoint, endPoint);
    isOutlier = false;
    density = 0.0;
}

CoarseSegment::CoarseSegment() {
    startPoint = QPointF(0.0, 0.0);
    endPoint = QPointF(0.0, 0.0);
    length = 0.0;
    isOutlier = false;
    maxAngle = 0.0;
    maxLength = 0.0;
    minLength= 0.0;
    maxLk = 0.0;
}

double CoarseSegment::calculateLength() {
    length = euclideanDistance(listFineSegments.first()->startPoint, listFineSegments.last()->endPoint);
    return length;
}

TrajectoryTRAOD::TrajectoryTRAOD()
{
    time = 0.0;
    score = 0.0;
    isOutlier = false;
}

void TrajectoryTRAOD::addPoint(QPointF point) {
    listPoints.append(point);
}

qint64 TrajectoryTRAOD::createSegments() {
    qint64 numOfSegment = 0;
    for (int i = 0; i < listPoints.size()-1; i++) {
        Segment *newSeg = new Segment(listPoints[i], listPoints[i+1]);
        listSegments.append(newSeg);
        numOfSegment++;
    }
    return numOfSegment;
}




