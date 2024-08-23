#ifndef GAMESCENE_H
#define GAMESCENE_H

#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <QObject>
#include "player.h"
#include "user.h"
#include "bug.h"
#include <QGraphicsView>
#include "deathscreen.h"
#include <QHBoxLayout>
#include <QSoundEffect>
#include <QLabel>

#pragma once



class GameScene : public QGraphicsScene
{
    Q_OBJECT
public:
    QGraphicsPixmapItem *pixMapItem;
    explicit GameScene(QObject *parent = nullptr);
    Player *player;
    User gameUser;
    int score;

    void takeDownGame();

public slots:
    void spawnBug();
    void startGameScene(User user, int levelNum);
    void onBugMissed();
    void onBugCaught();

private:
    int missedBugs;
    int caughtBugs;
    int maxMissedBugsAllowed;
    int level;
    int currentBugSpeed;
    int speedDoubled;
    QGraphicsView* gameView_;
    QWidget *mainWidget;
    deathScreen* gameOverScreen;
    QTimer *spawn_timer;
    QGraphicsView *mainView;
    QLabel *scoreLabel;

    int randomBugPosition();
    void setUserInfoStyle(QHBoxLayout *userInfoLayout);
    QSoundEffect *soundEffectMissed;
    QSoundEffect *soundEffectCaught;

signals:
    void gameOverSignal(int finalScore);
};

#endif // GAMESCENE_H
