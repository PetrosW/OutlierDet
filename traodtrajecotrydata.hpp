#ifndef TRAODTRAJECOTRYDATA_HPP
#define TRAODTRAJECOTRYDATA_HPP

#include <QObject>
#include <QPoint>
#include <QtMath>
#include <QFile>
#include <QDebug>

#include <traodtrajectory.hpp>

class TRAODTrajecotryData
{

public:
    TRAODTrajecotryData();

    // structure for holding data for detection
    QList<TrajectoryTRAOD *> listTrajectories;

    // clean lists so detection can start again with different parameters
    void cleanLists();

    // loading data (files) to trajectory structures
    // if you need load your own data create your own function
    // your new function must fill listPoints in TrajectoryIBoat object for every trajectory
    // and object of every trajectory add to listTrajectories
    bool readHurricaneFromFile(QString filePath);
    bool readPortoTaxiFromFile(QString filePath);


};

#endif // TRAODTRAJECOTRYDATA_HPP
