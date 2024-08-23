#ifndef DEATHSCREEN_H
#define DEATHSCREEN_H

#include <QDialog>
#include "user.h"

namespace Ui {
class deathScreen;
}

class deathScreen : public QDialog
{
    Q_OBJECT

public:
    explicit deathScreen(QWidget *parent = nullptr);
    ~deathScreen();
    void displayDeath(int currentScore, User user);

private slots:
    void on_replayButton_clicked();


private:
    Ui::deathScreen *ui;

signals:
    void replayClicked();

};

#endif // DEATHSCREEN_H
