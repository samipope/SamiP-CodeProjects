#include "player.h"

Player::Player(QObject *parent)
    : QObject{parent}
{
    setPixmap((QPixmap("://images/sami-hacker.png")).scaled(150,150));
    setPos(400, 365);
    setFlag(QGraphicsItem::ItemIsFocusable);

}

void Player::keyPressEvent(QKeyEvent *event){
    int inc = 40;

    int xPos = this->x();
    int yPos = this->y();

    // move right
    if(event->key() == Qt::Key_Right) {
        if(xPos < 800) {
            setPos(xPos + inc, yPos);
        }
    }

    // move left
    if(event->key() == Qt::Key_Left) {
        if(xPos > -20) {
            setPos(xPos - inc, yPos);
        }
    }
}
