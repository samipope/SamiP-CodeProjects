#include "login.h"
#include "ui_login.h"
#include "user.h"
#include "signup.h"
#include "gamescene.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>
#include <QDir>

#include <cctype> //for pass check

login::login(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::login)
    , signupWidget(new QWidget(this))  //create signup

{

    ui->setupUi(this);
    QPixmap logo(":/images/girlsJustWannaDebug.png");
    QPixmap scaledLogo = logo.scaled(ui->LogoLabel->width(), ui->LogoLabel->height(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    ui->LogoLabel->setPixmap(scaledLogo);


    QPixmap computerGirl(":/images/sami-hacker.png");
    QPixmap scaledComputerGirl = computerGirl.scaled(ui->imageLabel1->width(), ui->imageLabel1->height(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    ui->imageLabel1->setPixmap(scaledComputerGirl);


    QPixmap bug(":/images/glow-bug.png");
    QPixmap scaledBug = bug.scaled(ui->imageLabel2->width(), ui->imageLabel2->height(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    ui->imageLabel2->setPixmap(scaledBug);

    signupWidget->hide();  // Initially hide the signup widget


        this->setStyleSheet(

                        "QWidget {"
                        "background-color: rgba(0, 0, 0, 1);"
                        "}"

                        "QLineEdit {"
                        "background-color: rgba(0, 0, 0, 1);"
                        "color: #8D115F;"             // Text color is white
                        "border: 1px solid #8D115F;"
                        "padding: 2px;"
                        "font-family: 'Courier New';"
                        "}"

                        "QLabel {"
                        "color: white;"  // Change label text color to white for contrast
                        "font-family: 'Courier New';"
                       "font-weight: bold;"  // Make the font bold for better readability

                        "}"

                        "QPushButton#loginButton {"
                        "background-color: #8D115F;"  // Choose a suitable color
                        "color: white;"  // White text for contrast
                        "border: none;"  // Optional: remove border for a flat design
                        "padding: 5px;"  // Add padding for a bigger click area
                        "border-radius: 5px;"  // Rounded corners for the buttons
                        "font-family: 'Courier New';"
                        "font-weight: bold;"  // Make the font bold for better readability
                        "font-size: 35px;"
                        "}"

                        "QPushButton {"
                        "background-color: #EB5FB8;"  // Choose a suitable color
                        "color: #480831;"  // White text for contrast
                        "border: none;"  // Optional: remove border for a flat design
                        "padding: 5px;"  // Add padding for a bigger click area
                        "border-radius: 5px;"  // Rounded corners for the buttons
                        "font-family: 'Courier New';"
                        "font-weight: bold;"  // Make the font bold for better readability
                        "font-size: 20px;"
                        "}"

                        "QComboBox {"
                        "background-color: rgba(0, 0, 0, 1);" // Solid black background
                        "color: white;"                     // Bright pink text
                        "border: 1px solid #8D115F;"          // Dark pink border
                        "border-radius: 2px;"
                        "font-family: 'Courier New';"
                        "}"
                       );
}

login::~login()
{
    delete ui;
}


void login::on_playAsGuest_clicked()
{
    User guest;
    int theLevel = levelOn();
    emit loginSignal(guest.userName, theLevel);

}


void login::on_loginButton_clicked()
{

    QString qUsername = ui->usernameInput->text();
    QString qPassword = ui->passwordInput->text();
    std::string password = qPassword.toStdString(); //need this to use string func to chec


    if(qUsername.isEmpty() || qUsername.isNull()){
        QMessageBox::warning(this, "Username Empty", "username field cannot be empty");
    }


    //verify the password is correct
    if(User::checkPassword(qUsername,qPassword)){ //sucessful login
        //storedUserName = qUsername;
        int theLevel = levelOn();
        //qDebug() << "Level user selected: "+  QString::number(theLevel);

        //successful login, tell the controller
        emit loginSignal(qUsername, theLevel);
    }
    else{
        QMessageBox::warning(this, "Password and Username do not match", "Please try again or sign up");
        return;
    }
}


void login::on_signupButton_clicked()
{
    signup signup(this);
    if(signup.exec() == QDialog::Accepted){
        //if it is accepted, then we process the data
    }
}


int login::levelOn()
{

 QString level = ui->levelBox->currentText();
    //qDebug()<< level;

    if(level=="Level 2"){
     return 2; }
    else if(level=="Level 3"){
        return 3;}
    else return 1; //default is level one

}
