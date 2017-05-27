#include "traod.hpp"

TRAOD::TRAOD()
{
    MDLCostAdjustment = 10.0;
    numOfAllCoarseSegments = 0;
    numOfAllSegments = 0;

    p = 0.0;
    F = 0.0;
    D = 0.0;
    cl1 = new QList<Segment *>;
    cl2 = new QList<QList<Segment *> *>;
    wa = 0.0;
    wper = 0.0;
    wpar = 0.0;

    averageDensity = 0.0;
    sumDensity = 0.0;
}

bool TRAOD::loadTRAODHurricane(QString file) {
    return trajectoryData.readHurricaneFromFile(file);
}

bool TRAOD::loadTRAODPortoTaxi(QString file) {
    return trajectoryData.readPortoTaxiFromFile(file);
}

void TRAOD::runTRAOD() {
    numOfAllSegments = 0;
    foreach (auto traj, trajectoryData.listTrajectories) {
        numOfAllSegments += traj->createSegments();
    }
    createCoarseSegments();
    stdDev = standardDeviationAprox();
    densityForAllSegmenst();
    processOrigOpt();
}

void TRAOD::setParams(double p, double F, double D, double wa, double wper, double wpar) {
    this->p = p;
    this->F = F;
    this->D = D;
    this->wa = wa;
    this->wpar = wpar;
    this->wper = wper;
}

void TRAOD::processOrigOpt() {
    qint64 numOfTraj = trajectoryData.listTrajectories.size();
    qDebug() << "numOfTraj: " << numOfTraj;
    for (qint64 i = 0; i < trajectoryData.listTrajectories.size(); i++) {
        qint64 numOfOutlierSegments = 0;
        for (qint64 j = 0; j < trajectoryData.listTrajectories[i]->listCoarseSegments.size(); j++) {
            processSegmentOrigOpt(trajectoryData.listTrajectories[i]->listCoarseSegments[j], trajectoryData.listTrajectories[i]);
            QList<qint64> numOfCloseT = getNumOfCloseTOpt();
            // check all fine segments of their are anomalous
            for (qint64 k = 0; k < cl1->size(); k++) {
                double adj = computeAdjustingCoefficientOpt(cl1->at(k));
                if (qCeil(numOfCloseT[k] * adj) <= qCeil((1-p) * numOfTraj)) {
                    qDebug() << "Find Outlier Segment";
                    cl1->at(k)->isOutlier = true;
                    numOfOutlierSegments++;
                }
            }
        }
        if (isTrajectoryOutlier(&(trajectoryData.listTrajectories[i]->listSegments))) {
            qDebug() << "Find Outlier Trajectory";
            trajectoryData.listTrajectories[i]->isOutlier = true;
        }
        qDebug() << "traj out" << i << ":" << QDateTime::currentDateTime().time();
    }
}


void TRAOD::processSegmentOrigOpt(CoarseSegment *coarseSegment, TrajectoryTRAOD *trajIn) {
    cl1->clear();
    cl2->clear();
    foreach (const TrajectoryTRAOD *traj, trajectoryData.listTrajectories) {
        if (traj == trajIn) continue;
        QList<Segment *> *clTraj2 = new  QList<Segment *>();
        foreach (CoarseSegment *coarseSeg, traj->listCoarseSegments) {
            QPair<double, double> bounds = getBounds(coarseSegment, coarseSeg);
            if (bounds.first > D) {
                continue;
            }
            if (bounds.second <= D) {
                foreach (Segment *fineSeg, coarseSegment->listFineSegments) {
                    if (!cl1->contains(fineSeg)) {
                        cl1->append(fineSeg);
                    }
                }
                foreach (Segment *fineSeg, coarseSeg->listFineSegments) {
                    if (!clTraj2->contains(fineSeg)) {
                        clTraj2->append(fineSeg);
                    }
                }
            }
            else {
                foreach (Segment *fineSeg1, coarseSegment->listFineSegments) {
                    foreach (Segment *fineSeg2, coarseSeg->listFineSegments) {
                        if (disAllFin(fineSeg1, fineSeg2) < D) {
                            if (!cl1->contains(fineSeg1)) {
                                cl1->append(fineSeg1);
                            }
                            if (!clTraj2->contains(fineSeg2)) {
                                clTraj2->append(fineSeg2);
                            }
                        }

                    }
                }
            }
        }
        if (!clTraj2->isEmpty()) {
            cl2->append(clTraj2);
        }
    }
}

QPair<double, double> TRAOD::getBounds(CoarseSegment *coarseSeg1, CoarseSegment *coarseSeg2) {
    double k1 = 0.0;
    double k2 = 0.0;
    double lk1 = 0.0;
    double lk2 = 0.0;
    QPair<double, double> lp(0.0, 0.0);
    double angle = 0.0;
    if (coarseSeg1->length > coarseSeg2->length) {
        QPointF p1 = perpenPoint(coarseSeg1->startPoint, coarseSeg1->endPoint, coarseSeg2->startPoint, &k1);
        QPointF p2 = perpenPoint(coarseSeg1->startPoint, coarseSeg1->endPoint, coarseSeg2->endPoint, &k2);
        lp = disPar(k1, k2, coarseSeg1->startPoint, coarseSeg1->endPoint, p1, p2);
        lk1 = euclideanDistance(p1, coarseSeg2->startPoint);
        lk2 = euclideanDistance(p2, coarseSeg2->endPoint);
        double L = disAngle(coarseSeg2->length, coarseSeg1->startPoint, coarseSeg1->endPoint,  coarseSeg2->startPoint, coarseSeg2->endPoint, &angle);
    }
    else {

        QPointF p1 = perpenPoint(coarseSeg2->startPoint, coarseSeg2->endPoint, coarseSeg1->startPoint, &k1);
        QPointF p2 = perpenPoint(coarseSeg2->startPoint, coarseSeg2->endPoint, coarseSeg1->endPoint, &k2);
        lp = disPar(k1, k2, coarseSeg2->startPoint, coarseSeg2->endPoint, p1, p2);
        lk1 = euclideanDistance(p1, coarseSeg1->startPoint);
        lk2 = euclideanDistance(p2, coarseSeg1->endPoint);
        double L = disAngle(coarseSeg1->length, coarseSeg2->startPoint, coarseSeg2->endPoint,  coarseSeg1->startPoint, coarseSeg1->endPoint, &angle);
    }
    double lowerBoundPer = qMin(lk1, lk2) - (coarseSeg1->maxLk + coarseSeg2->maxLk);
    double upperBoundPer = qMax(lk1, lk2) + (coarseSeg1->maxLk + coarseSeg2->maxLk);

    double dPar = qMin(lp.first, lp.second);

    double lowerBoundPar = 0.0;
    double upperBoundPar = 0.0;

    if (qFabs(k1) < 1.0  && qFabs(k2) < 1.0) {  // enclose
        upperBoundPar = qMax(coarseSeg1->length, coarseSeg2->length);
    }
    else if (qFabs(k1) < 1.0  ||  qFabs(k2) < 1.0) { // overlap
        upperBoundPar = coarseSeg1->length + coarseSeg2->length - dPar;
    }
    else { // disjoint
        lowerBoundPar = dPar;
        upperBoundPar = coarseSeg1->length + coarseSeg2->length + dPar;
    }

    double minAngle = (angle - coarseSeg1->maxAngle - coarseSeg2->maxAngle);
    double maxAngle = (angle + coarseSeg1->maxAngle + coarseSeg2->maxAngle);

    if (minAngle < 0.0) minAngle = 0.0;
    if (maxAngle > M_PI_2) maxAngle = M_PI_2;

    double lowerBoundAngle = qMin(coarseSeg1->minLength, coarseSeg2->minLength) * qSin(minAngle);
    double upperBoundAngle = qMin(coarseSeg1->maxLength, coarseSeg2->maxLength) * qSin(maxAngle);

    double lowerBound = weightDis(lowerBoundAngle, lowerBoundPer, lowerBoundPar);
    double upperBound = weightDis(upperBoundAngle, upperBoundPer, upperBoundPar);
    return QPair<double, double>(lowerBound, upperBound);

}

QList<qint64> TRAOD::getNumOfCloseTOpt() {
    QList<qint64> numOfCloseTL;
    foreach (Segment *seg1, *cl1) {
        qint64 numOfCloseT = 0;
        foreach (auto traj, *cl2) {
            double lengthOfSegments = 0.0;
            foreach (Segment *seg2, *traj) {
                if (seg1 == seg2) continue;
                if (disAllFin(seg1, seg2) < D) {
                    lengthOfSegments += seg2->length;
                    if (seg1->length < lengthOfSegments) {
                        numOfCloseT++;
                        break;
                    }
                }
            }
        }
        numOfCloseTL.append(numOfCloseT);
    }
    return numOfCloseTL;
}

double TRAOD::computeAdjustingCoefficientOpt(Segment *seg) {
    if (seg->density == 0.0) {
        return 0.0;
    }
    return averageDensity / seg->density;
}

bool TRAOD::isTrajectoryOutlier(QList<Segment *> *listSegments) {
    double normalLength = 0.0;
    double anomalyLength = 0.0;
    foreach (Segment *seg, *listSegments) {
        if (seg->isOutlier) {
            anomalyLength += seg->length;
        }
        normalLength += seg->length;
    }
    if (normalLength == 0.0)
        return true;
    if ((anomalyLength / normalLength) >= F)
        return true;
    return false;
}

void TRAOD::createCoarseSegments() {
    foreach (TrajectoryTRAOD *traj, trajectoryData.listTrajectories) {
        qint64 startIndex = 0;
        qint64 length = 1;
        traj->listCoarsePoints.append(traj->listPoints[startIndex]);
        CoarseSegment *newCoarseSeg = new CoarseSegment();
        newCoarseSeg->listFineSegments.append(traj->listSegments[startIndex]);
        newCoarseSeg->startPoint = traj->listSegments[startIndex]->startPoint;
        while (startIndex + length < traj->listSegments.size()) {
            qint64 currentIndex = startIndex + length;
            QVector<double> *coarseSegParam = new QVector<double>(4, 0.0);
            double costNoPar = MDLnopar(traj->listSegments[startIndex], traj->listSegments[currentIndex]);
            double costPar = MDLpar(startIndex, currentIndex, &(traj->listSegments), *coarseSegParam);
            if (costPar > (costNoPar + MDLCostAdjustment)) {
                // if consPar is higher last is add to coarse segment
                // and new coarse segment object is created
                traj->listCoarsePoints.append(traj->listPoints[currentIndex - 1]);
                newCoarseSeg->calculateLength();
                traj->listCoarseSegments.append(newCoarseSeg);
                numOfAllCoarseSegments++;
                newCoarseSeg = new CoarseSegment();
                startIndex = currentIndex;
                length = 1;
                // and start point and start fine segment to new coarse segment
                newCoarseSeg->listFineSegments.append(traj->listSegments[startIndex]);
                newCoarseSeg->startPoint = traj->listSegments[startIndex]->startPoint;
            }
            else {
                // if consPar is lower fine segment is add to coarse segment
                length += 1;
                newCoarseSeg->listFineSegments.append(traj->listSegments[currentIndex]);
                newCoarseSeg->endPoint = traj->listSegments[currentIndex]->endPoint;
                newCoarseSeg->maxLength = coarseSegParam->at(0);
                newCoarseSeg->minLength = coarseSegParam->at(1);
                newCoarseSeg->maxAngle = coarseSegParam->at(2);
                newCoarseSeg->maxLk = coarseSegParam->at(3);
            }
        }
        traj->listCoarsePoints.append(traj->listPoints.last());
        newCoarseSeg->listFineSegments.append(traj->listSegments.last());
        newCoarseSeg->endPoint = traj->listSegments.last()->endPoint;
        newCoarseSeg->calculateLength();
        traj->listCoarseSegments.append(newCoarseSeg);
        numOfAllCoarseSegments++;
    }
}

double TRAOD::MDLnopar(Segment *L1, Segment *L2) {
    double ed = euclideanDistance(L1->startPoint, L2->endPoint);
    // because of log2 (log2(1.0) == 0)
    if (ed <= 0.0) {
        ed = 1.0;
    }
    double MDL = log2(ed);
    return qFabs(MDL);
}


//newCoarseSeg->maxLength = coarseSegParam[0];
//newCoarseSeg->minLength = coarseSegParam[1];
//newCoarseSeg->maxAngle = coarseSegParam[2];
//newCoarseSeg->maxPerpen = coarseSegParam[3];
double TRAOD::MDLpar(qint64 startIndex, qint64 endIndex, QList<Segment *>* listSegments, QVector<double> &coarseSegParam) {
    double MDLcon = MDLnopar(listSegments->at(startIndex), listSegments->at(endIndex));
    QPointF startPartSeg = listSegments->at(startIndex)->startPoint;
    QPointF endPartSeg = listSegments->at(endIndex)->endPoint;
    double partSegLength = euclideanDistance(startPartSeg, endPartSeg);
    coarseSegParam[1] = partSegLength;
    for (qint64 i = startIndex; i < endIndex; ++i) {
        QPointF startPoint = listSegments->at(i)->startPoint;
        QPointF endPoint = listSegments->at(i)->endPoint;
        double k1 = 0.0;
        double k2 = 0.0;
        QPointF p1 = perpenPoint(startPartSeg, endPartSeg, startPoint, &k1);
        QPointF p2 = perpenPoint(startPartSeg, endPartSeg, endPoint, &k2);
        double lk1 = euclideanDistance(p1, startPoint);
        double lk2 = euclideanDistance(p2, endPoint);
        double maxlk = qMax(lk1, lk2);
        double dPer = disPerpen(lk1, lk2);
        double segLength = listSegments->at(i)->length;
        double angle = 0.0;
        double dAngle = disAngle(segLength, startPartSeg, endPartSeg, startPoint, endPoint, &angle);
        if (dPer < 1.0) { dPer = 1.0; }
        if (dAngle < 1.0) { dAngle = 1.0; }

        if (coarseSegParam[0] < segLength) coarseSegParam[0] = segLength;
        if (coarseSegParam[1] < segLength) coarseSegParam[1] = segLength;
        if (coarseSegParam[2] < angle) coarseSegParam[2] = angle;
        if (coarseSegParam[3] < maxlk) coarseSegParam[3] = maxlk;
        MDLcon += qFabs(log2(dPer)) + qFabs(log2(dAngle));
    }
    return MDLcon;
}

double TRAOD::weightDis(double disAngle, double disPer, double disPar) {
    return disAngle * wa + disPer * wper + disPar * wpar;
}

double TRAOD::disAllFin(Segment *L1, Segment *L2) {

    double lp1 = 0.0;
    double lp2 = 0.0;
    double lk1 = 0.0;
    double lk2 = 0.0;
    double da = 0.0;
    double angle = 0.0;
    if (L1->length == 0.0 || L2->length == 0.0) {
        lk1 = euclideanDistance(L1->startPoint, L2->startPoint);
        lk2 = euclideanDistance(L1->endPoint, L2->endPoint);
    }
    else {
        if (L1->length > L2->length) {
            double k1 = 0.0;
            double k2 = 0.0;
            QPointF p1 = perpenPoint(L1->startPoint, L1->endPoint, L2->startPoint, &k1);
            QPointF p2 = perpenPoint(L1->startPoint, L1->endPoint, L2->endPoint, &k2);
            QPair<double, double> lp = disPar(k1, k2, L1->startPoint, L1->endPoint, p1, p2);
            lp1 = lp.first;
            lp2 = lp.second;
            lk1 = euclideanDistance(p1, L2->startPoint);
            lk2 = euclideanDistance(p2, L2->endPoint);
            da = disAngle(L2->length, L1->startPoint, L1->endPoint, L2->startPoint, L2->endPoint, &angle);
        }
        else {
            double k1 = 0.0;
            double k2 = 0.0;
            QPointF p1 = perpenPoint(L2->startPoint, L2->endPoint, L1->startPoint, &k1);
            QPointF p2 = perpenPoint(L2->startPoint, L2->endPoint, L1->endPoint, &k2);
            QPair<double, double> lp = disPar(k1, k2, L2->startPoint, L2->endPoint, p1, p2);
            lp1 = lp.first;
            lp2 = lp.second;
            lk1 = euclideanDistance(p1, L1->startPoint);
            lk2 = euclideanDistance(p2, L1->endPoint);
            da = disAngle(L1->length, L2->startPoint, L2->endPoint, L1->startPoint, L1->endPoint, &angle);
        }
    }
    double dper = disPerpen(lk1, lk2);
    double dpar = qMin(lp1, lp2);
    return weightDis(da, dper, dpar);
}


QPointF TRAOD::perpenPoint(QPointF LS, QPointF LE, QPointF LP, double *k) {
    // create vectors
    QVector2D v1(LE.x() - LS.x(), LE.y() - LS.y());
    QVector2D v2(LP.x() - LS.x(), LP.y() - LS.y());
    if (v1[0] == 0.0 && v1[1] == 0.0) {
        *k = 0.0;
        return QPointF(0.0, 0.0);
    }
    // inner product / length of vectors
    *k = (v1[0] * v2[0] + v1[1] * v2[1]) / (qPow(v1[0],2) + qPow(v1[1],2));
    QPointF point(LS.x() + *k * v1[0], LS.y() + *k * v1[1]);
    return point;
}

double TRAOD::disAngle(double L2, QPointF L1S, QPointF L1E, QPointF L2S, QPointF L2E, double *angle) {
    // create vectors
    QVector2D v1(L1E.x() - L1S.x(), L1E.y() - L1S.y());
    QVector2D v2(L2E.x() - L2S.x(), L2E.y() - L2S.y());
    // inner product of the two vectors
    double d = (v1[0] * v2[0]) + (v1[1] * v2[1]);
    // length of vectors
    double l = qSqrt(qPow(v1[0],2) + qPow(v1[1],2)) * qSqrt(qPow(v2[0],2) + qPow(v2[1],2));
    if (l == 0.0) {
        *angle = 0.0;
        return 0.0;
    }
    double cosTheta = d / l;
    // compensate the computation error (e.g., 1.00001)
    if (cosTheta > 1.0) cosTheta = 1.0;
    if (cosTheta < -1.0) cosTheta = -1.0;
    double sinTheta = qSqrt(1 - qPow(cosTheta, 2));
    *angle = qAsin(sinTheta);
    return (L2 * sinTheta);
}

double TRAOD::disPerpen(double lk1, double lk2) {
    if (lk1 == 0.0 && lk2 == 0.0) {
        return 0.0;
    }
    return (qPow(lk1, 2) + qPow(lk2, 2)) / (lk1 + lk2);
}

QPair<double, double> TRAOD::disPar(double k1, double k2, QPointF S, QPointF E, QPointF p1, QPointF p2) {
    double lp1 = 0.0;
    double lp2 = 0.0;
    if (k1 < 0.5)
        lp1 = euclideanDistance(S, p1);
    else
        lp1 = euclideanDistance(E, p1);
    if (k2 < 0.5)
        lp2 = euclideanDistance(S, p2);
    else
        lp2 = euclideanDistance(E, p2);
    return QPair<double, double>(lp1, lp2);

}

double TRAOD::standardDeviation() {
    std::list<double> listDist;
    qint64 k = 1;
    qDebug() << "numTraj" << trajectoryData.listTrajectories.size();
    foreach (TrajectoryTRAOD *traj1, trajectoryData.listTrajectories) {
        foreach (Segment *seg1, traj1->listSegments) {
            for (qint64 i = k; i < trajectoryData.listTrajectories.size(); i++) {
                foreach (Segment *seg2, trajectoryData.listTrajectories.at(i)->listSegments) {
                    listDist.push_back(disAllFin(seg1, seg2));
                }
            }
        }
        qDebug() << listDist.size();
        qDebug() << "std: " << k << ": " << QDateTime::currentDateTime().time();
        k++;
    }
    qint64 listDistSize = listDist.size();
    double mean = 0.0;
    qDebug() << "mean";
    foreach (double dist, listDist) {
        mean += (dist / listDistSize);
    }
    qDebug() << mean;
    std::list<double> listDeviation;
    qDebug() << "listDeviation";
    foreach (double dist, listDist) {
        listDeviation.push_back(qPow(dist - mean,2));
    }
    double stdDev = 0.0;
    qDebug() << "stdDev";
    qint64 listDeviationSize = listDeviation.size();
    foreach (double dev, listDeviation) {
        stdDev += (dev / listDeviationSize);
    }
    qDebug() << "stdDev" << stdDev;
    qDebug() << "qSqrt(stdDev)" << qSqrt(stdDev);
    return qSqrt(stdDev);
}

double TRAOD::standardDeviationAprox() {
    QList<qint64> listRandomIndex;
    qint64 i = 0;
    while(i < trajectoryData.listTrajectories.size()/10) {
        qint64 r = qrand() % (trajectoryData.listTrajectories.size());
        if (!listRandomIndex.contains(r)) {
            listRandomIndex.append(r);
            i++;
        }
    }

    qint64 k = 1;
    std::list<double> listDist;
    qDebug() << "numTraj" << trajectoryData.listTrajectories.size();
    foreach (qint64 r, listRandomIndex) {
        foreach (Segment *seg1, trajectoryData.listTrajectories[r]->listSegments) {
            for (qint64 i = 0; i < trajectoryData.listTrajectories.size(); i++) {
                if (r == i) continue; // not between same trajectory
                foreach (Segment *seg2, trajectoryData.listTrajectories.at(i)->listSegments) {
                    listDist.push_back(disAllFin(seg1, seg2));
                }
            }
        }
        qDebug() << listDist.size();
        qDebug() << "std: " << k << ": " << QDateTime::currentDateTime().time();
        k++;
    }
    qint64 listDistSize = listDist.size();
    double mean = 0.0;
    qDebug() << "mean";
    foreach (double dist, listDist) {
        mean += (dist / listDistSize);
    }
    qDebug() << mean;
    std::list<double> listDeviation;
    qDebug() << "listDeviation";
    foreach (double dist, listDist) {
        listDeviation.push_back(qPow(dist - mean,2));
    }
    double stdDev = 0.0;
    qDebug() << "stdDev";
    qint64 listDeviationSize = listDeviation.size();
    foreach (double dev, listDeviation) {
        stdDev += (dev / listDeviationSize);
    }
    qDebug() << "stdDev" << stdDev;
    qDebug() << "qSqrt(stdDev)" << qSqrt(stdDev);
    return qSqrt(stdDev);
}

void TRAOD::densityForAllSegmenst() {
    qint64 k = 0;
    sumDensity = 0;
    foreach (TrajectoryTRAOD *traj1, trajectoryData.listTrajectories) {
        foreach (Segment *seg1, traj1->listSegments) {
            qint64 numOfCloseS = 0;
            foreach (TrajectoryTRAOD *traj2, trajectoryData.listTrajectories) {
                if (traj1 == traj2) continue; // not between same trajectory
                foreach (Segment *seg2, traj2->listSegments) {
                    if (disAllFin(seg1, seg2) < stdDev) {
                        numOfCloseS++;
                    }
                }
            }
            seg1->density = numOfCloseS;
            sumDensity += numOfCloseS;
        }
        qDebug() << "densi " << k << ": " << QDateTime::currentDateTime().time();
        k++;
    }
    qDebug() << "sumDensity / numOfAllSegments: " << sumDensity / numOfAllSegments;
    averageDensity = sumDensity / numOfAllSegments;
}

double euclideanDistance(QPointF startPoint, QPointF endPoint) {
    return qSqrt( qPow(endPoint.x()-startPoint.x(), 2) + qPow(endPoint.y()-startPoint.y(), 2) );
}
