#include "bug.h"
#include "QtWidgets/qgraphicsscene.h"
#include <QSoundEffect>

int Bug::missedBugs = 0;

Bug::Bug(QObject *parent, int speed)
    : QObject{parent}, speed(speed)
{
    // set image
    setPixmap((QPixmap("://images/glow-bug.png")).scaled(30,30));
    //light pink: d9adc8
    //med pink: e56db7
    //dark pink: ea249f

    QTimer *move_timer = new QTimer(this);
    connect(move_timer, &QTimer::timeout, this, &Bug::moveBug);
    move_timer->start(speed);
}

void Bug::moveBug() {
    int inc = 2;
    int yPos = this->y();
    int xPos = this->x();

    // check if out of bounds, and delete if needed
    if(yPos > 510) {
        emit bugMissed();
        scene()->removeItem(this);
        delete this;
        return;
    }

    // check collision with player, and delete if needed
    if(yPos > 365) {
        QList<QGraphicsItem*> itemsThatCollide = scene()->collidingItems(this);
        for (QGraphicsItem* item : itemsThatCollide)
        {
            if (item->type() == Player::Type)
            {
                // Emit a signal to GameScene
                emit bugCaught();
                scene()->removeItem(this);
                delete this;
                return;
            }
        }
    }


    // move down
    setPos(xPos, yPos + inc);
    \
}


