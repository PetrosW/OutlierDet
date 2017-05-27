#include "iboattrajectory.hpp"

TrajectoryIBoat::TrajectoryIBoat(qint64 ID) : ID(ID) {
    time = 0.0;
    score.append(0.0);
}


CellIBoat::CellIBoat(QPoint cell) : cell(cell) {

}

void CellIBoat::addTraj(qint64 pos, TrajectoryIBoat *traj) {
    if (listTrajectory.contains(pos)) {
        listTrajectory[pos].append(traj);
    }
    else {
        QList<TrajectoryIBoat *> newList;
        newList.append(traj);
        listTrajectory.insert(pos, newList);
    }
}
