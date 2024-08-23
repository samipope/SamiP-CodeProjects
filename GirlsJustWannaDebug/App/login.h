#ifndef LOGIN_H
#define LOGIN_H


#include <QMainWindow>

namespace Ui {
class login;
}

class login : public QMainWindow
{
    Q_OBJECT

public:
    explicit login(QWidget *parent = nullptr);
    ~login();
    int levelOn();


signals:
    void loginSignal(const QString &data, int levelNum);


private slots:


    void on_playAsGuest_clicked();

    void on_loginButton_clicked();

    void on_signupButton_clicked();



private:
    Ui::login *ui;
    QWidget *signupWidget;  // Widget to hold the signup form
    // QWidget *deathWidget;
    //QString storedUserName;

};

#endif // LOGIN_H
