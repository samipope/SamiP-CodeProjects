#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "login.h"
#include "user.h"
#include "gamescene.h"
#include "deathscreen.h"

class Controller : public QGraphicsScene
{
    int level;
public:
    Controller();
    void setUpLogin();
    login loginScreen;
    User* player;
    GameScene* game;
    deathScreen* gameOverScreen;


public slots:
    void handleLogin(const QString &data, int levelNum);
    void handleGameOver(int value);
    void handleReplay();
};

#endif // CONTROLLER_H

