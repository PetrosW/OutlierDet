#include "traodtrajecotrydata.hpp"

TRAODTrajecotryData::TRAODTrajecotryData() {
//    listTrajectories = new QList<Trajectory *>();
}

bool TRAODTrajecotryData::readHurricaneFromFile(QString filePath) {
    cleanLists();
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QFile::Text)) {
        qDebug() << "error load file";
        return false;
    }
    file.readLine();
    file.readLine();
    while(!file.atEnd()) {
        TrajectoryTRAOD *newTraj = new TrajectoryTRAOD();
        QString line = file.readLine();
        QStringList data = line.split(" ");
        for (qint64 i = 2; i < data.size()-1; i+=2) {
            newTraj->addPoint(QPointF(data[i].toDouble(), data[i+1].toDouble()));
        }
        listTrajectories.append(newTraj);
    }
    file.close();
    return true;
}

bool TRAODTrajecotryData::readPortoTaxiFromFile(QString filePath) {
    cleanLists();
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QFile::Text)) {
        qDebug() << "error load file";
        return false;
    }
    while(!file.atEnd()) {
        TrajectoryTRAOD *newTraj = new TrajectoryTRAOD();
        QString line = file.readLine();
        QStringList data = line.split(";");
        newTraj->time = data[0].toDouble();
        for (qint64 i = 1; i < data.size(); i++) {
            QStringList data2 = data[i].split(",");
            double x = (data2[0].toDouble()-(-8.8))/((-8.5)-(-8.8)) * 1000;
            double y = (data2[1].toDouble()-41.09)/(41.27-41.09) * 1000;
            newTraj->addPoint(QPointF(x, y));
        }
        listTrajectories.append(newTraj);
    }
    file.close();
    return true;
}

void TRAODTrajecotryData::cleanLists() {
    foreach (TrajectoryTRAOD *traj, listTrajectories) {
        qDeleteAll(traj->listCoarseSegments.begin(), traj->listCoarseSegments.end());
        traj->listCoarseSegments.clear();
        qDeleteAll(traj->listSegments.begin(), traj->listSegments.end());
        traj->listSegments.clear();
    }
    qDeleteAll(listTrajectories.begin(), listTrajectories.end());
    listTrajectories.clear();
}
