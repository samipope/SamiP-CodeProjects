#include <QApplication>
#include <QGraphicsView>
#include <QGraphicsPixmapItem>
#include "login.h"
#include "user.h"
#include "gamescene.h"
#include "controller.h"




int main(int argc, char **argv) {
    QApplication app (argc, argv);

    // Load the icon from a file
    QIcon icon("://images/sami-hacker.png");

    // Set the application icon
    app.setWindowIcon(icon);

    //make a controller and have it set up the login screen
    Controller controller;
    controller.setUpLogin();

    return app.exec();
}



