#ifndef BUG_H
#define BUG_H

#include <QObject>
#include <QGraphicsPixmapItem>
#include <QTimer>
#include "player.h"

class Bug : public QObject, public QGraphicsPixmapItem
{
    Q_OBJECT
private:
    int speed;
public:
    explicit Bug(QObject *parent = nullptr, int speed = 100);
    static int missedBugs;

public slots:
    void moveBug();

signals:
    void bugMissed();
    void bugCaught();



};

#endif // BUG_H
