#include "mainwindow.h"
#include <QApplication>
#include <QTime>
#include <QString>
#include <QDesktopWidget>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
        app.setApplicationDisplayName("TSP-Solver v1.0");
        app.setApplicationName("Traveler Salesman Problem Solver v1.0");
        app.setApplicationVersion("1.0");
        app.setWindowIcon(QIcon("icon.ico"));

    qsrand(static_cast<unsigned>(QTime(0,0,0).secsTo(QTime::currentTime())));
    QFont tfont("YouYuan",13);
    app.setFont(tfont);
    MainWindow w;

    w.show();

    return app.exec();
}
