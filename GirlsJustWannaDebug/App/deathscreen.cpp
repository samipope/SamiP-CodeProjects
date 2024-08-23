#include "deathscreen.h"
#include "user.h"
#include "ui_deathscreen.h"

deathScreen::deathScreen(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::deathScreen)
{
    ui->setupUi(this);





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
        "color: #8D115F;"  // Change label text color to white for contrast
        "font-family: 'Courier New';"
        "font-weight: bold;"  // Make the font bold for better readability

        "}"

        "QPushButton {"
        "background-color: #EB5FB8;"  // Choose a suitable color
        "color: #480831;"  // White text for contrast
        "border: none;"  // Optional: remove border for a flat design
        "padding: 5px;"  // Add padding for a bigger click area
        "border-radius: 5px;"  // Rounded corners for the buttons
        "font-family: 'Courier New';"
        "font-weight: bold;"  // Make the font bold for better readability

        "}"

        "QComboBox {"
        "background-color: rgba(0, 0, 0, 1);" // Solid black background
        "color: white;"
        "border: 1px solid #8D115F;"          // Dark pink border
        "padding: 2px;"
        "border-radius: 2px;"
        "font-family: 'Courier New';"
        "font-weight: bold;"  // Make the font bold for better readability

        "}"



        );
}



deathScreen::~deathScreen()
{
    delete ui;
}


//IF REPLAY BUTTON IS PRESSED
void deathScreen::on_replayButton_clicked()
{
    emit replayClicked();
}



//maybe this needs to take the user instead/in addition to?
void deathScreen::displayDeath(int currentScore, User user){

    //qDebug() << "reached deathScreen class";


    if(currentScore>=100){
        QPixmap winImage(":/images/wonImage.png");
        QPixmap scaledWinImage = winImage.scaled(ui->deathImageLabel->width(), ui->deathImageLabel->height(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
        ui->deathImageLabel->setPixmap(scaledWinImage);

    }
    else{
        QPixmap deathImage(":/images/deathImage.png");
        QPixmap scaledDeathImage = deathImage.scaled(ui->deathImageLabel->width(), ui->deathImageLabel->height(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
        ui->deathImageLabel->setPixmap(scaledDeathImage);
    }



    if(QString::compare(user.userName, "guest", Qt::CaseInsensitive) == 0){
        ui->yourHighScore->setText("Set up an account to store your high score!");

    }
    else{ //set and presents high score if it is a valid user
        ui->yourHighScore->setText("Your High score: " + QString::number(user.getUserHighScore()) );
    }


    //print the final score
    QString score = QString::number(currentScore);
    ui->yourScore->setText("Your score: " + score);

    //print the global high score
    ui->globalHighScore->setText("Global High Score: " + QString::number(user.getGlobalHighScore()));

    //show the death screen
    this->show();
}

