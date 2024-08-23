#ifndef PLAYER_H
#define PLAYER_H

#include <QObject>
#include <QGraphicsPixmapItem>
#include <QKeyEvent>

class Player : public QObject, public QGraphicsPixmapItem
{
    Q_OBJECT
public:
    explicit Player(QObject *parent = nullptr);

    void keyPressEvent(QKeyEvent *event);

signals:
};

#endif // PLAYER_H
