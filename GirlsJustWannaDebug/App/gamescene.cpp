#include "gamescene.h"
#include <QGraphicsView>
#include <QLabel>
#include <QHBoxLayout>
#include <QApplication>


GameScene::GameScene(QObject *parent)
    : QGraphicsScene{parent}
{

    //set background
    setSceneRect(0,0,908,510);
    setBackgroundBrush(QBrush(QImage("://images/binaryBackground.jpg").scaledToHeight(512) .scaledToWidth(910)));

    // add player (bucket)
    player = new Player();
    player->setFocus();
    addItem(player);

    // set bug counter and score
    missedBugs = 0;
    caughtBugs = 0;
    maxMissedBugsAllowed = 5;
    score = 0;

    // set up spawn timer
    level = 1;
    currentBugSpeed = 200;
    speedDoubled = 0;
    spawn_timer = new QTimer(this);
    connect(spawn_timer, &QTimer::timeout, this, &GameScene::spawnBug);
    spawn_timer->start(2000);



    //sound effects
    soundEffectMissed = new QSoundEffect(this);
    soundEffectMissed->setSource(QUrl::fromLocalFile("://sounds/missedBug.wav"));
    soundEffectMissed->setVolume(0.25);

    soundEffectCaught = new QSoundEffect(this);
    soundEffectCaught->setSource(QUrl::fromLocalFile("://sounds/caughtBug.wav"));
    soundEffectCaught->setVolume(0.25);
}


int GameScene::randomBugPosition() {
    // get random x coord: 0 < x < 900
    int random_number = arc4random() % 901;
    return random_number;
}

void GameScene::spawnBug() {
    Bug *bug = new Bug(nullptr, currentBugSpeed);
    bug->setPos(randomBugPosition(), 0);
    addItem(bug);
    connect(bug, &Bug::bugMissed, this, &GameScene::onBugMissed);
    connect(bug, &Bug::bugCaught, this, &GameScene::onBugCaught);
}

void GameScene::startGameScene(User user, int levelNum){
    //store the user data
    gameUser = user;
    //set the level and bug speed
    level = levelNum;
    currentBugSpeed /= level;

    // Create the main widget to hold the QGraphicsView and user info
    mainWidget = new QWidget();
    QVBoxLayout *mainLayout = new QVBoxLayout(mainWidget);

    // Create a label for the user name
    QLabel *userNameLabel = new QLabel(user.userName);

    // Create a label for the profile picture
    QLabel *profilePicLabel = new QLabel;
    QPixmap profilePic(user.profilePicPath);

    // Scale the profile picture to a specific number of pixels
    int targetWidth = 50; // Set the target width in pixels
    int targetHeight = 50; // Set the target height in pixels
    QPixmap scaledProfilePic = profilePic.scaled(targetWidth, targetHeight, Qt::KeepAspectRatio);
    profilePicLabel->setPixmap(scaledProfilePic);

    // Create a label for the current date
    QLabel *dateLabel = new QLabel;
    QString currentDate = QDateTime::currentDateTime().toString("MM-dd-yyyy"); // Get the current date in the desired format
    dateLabel->setText("\tDate: " + currentDate);

    // Create a label for the player's high score
    scoreLabel = new QLabel;
    scoreLabel->setText("\tScore: " + QString::number(score));


    // Add the user name and profile picture labels to a horizontal layout
    QHBoxLayout *userInfoLayout = new QHBoxLayout;
    userInfoLayout->addWidget(profilePicLabel);
    userInfoLayout->addWidget(userNameLabel);
    userInfoLayout->addWidget(dateLabel);
    userInfoLayout->addWidget(scoreLabel);
    setUserInfoStyle(userInfoLayout);

    // Set the alignment of the user info layout to top left
    userInfoLayout->setAlignment(Qt::AlignTop | Qt::AlignLeft);

    // Add the user info layout to the main layout
    mainLayout->addLayout(userInfoLayout);

    // Create the QGraphicsView for the game scene
    mainView = new QGraphicsView();
    mainView->setFixedSize(910, 512);
    mainView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    mainView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    // Set the scene to the view
    mainView->setScene(this);

    // Add the QGraphicsView to the main layout
    mainLayout->addWidget(mainView);

    // Show the main widget
    // Get the geometry of the primary screen
    QRect screenGeometry = QGuiApplication::primaryScreen()->geometry();

    // Calculate the center position of the widget
    int x = (screenGeometry.width() - 910) / 2;
    int y = (screenGeometry.height() - 700) / 2;
    mainWidget->move(x,y);
    mainWidget->show();

    // Set focus to the bucket object
    player->setFocus();
}


void GameScene::onBugMissed()
{
    missedBugs++;
    if (missedBugs >= maxMissedBugsAllowed) {
        //update the scores
        gameUser.updateScores(score);
        gameUser.save();
        //tell the controller there was a game over
        emit gameOverSignal(score);
    } else{
        soundEffectMissed->play();
    }

}

void GameScene::onBugCaught()
{
    score+=5;
    caughtBugs++;
    soundEffectCaught->play();
    scoreLabel->setText("\tScore: " + QString::number(score));

    if(score >= 150) {
        gameUser.updateScores(score);
        gameUser.save();
        emit gameOverSignal(score);
    }


    if(caughtBugs % 5 == 0 && speedDoubled < 5) {
        speedDoubled++;
        currentBugSpeed = currentBugSpeed / 2;
    }

}


void GameScene::setUserInfoStyle(QHBoxLayout *userInfoLayout){
    // Iterate through each widget in the layout and set styles for QLabel
    for (int i = 0; i < userInfoLayout->count(); ++i) {
        QWidget* widget = userInfoLayout->itemAt(i)->widget();
        if (widget) {
            QLabel* label = qobject_cast<QLabel*>(widget);
            if (label) {
                label->setStyleSheet("color: #EB5FB8; font-family: 'Courier New'; font-weight: bold; font-size: 16px;");
            }
        }
    }
}

// handles the closing of the game window
void GameScene::takeDownGame(){
    //qDebug() <<"called takeDownGame";

    //stop spawning bugs
    if(this->spawn_timer){
        //qDebug() << "stopping spawn timer";
        spawn_timer->stop();
    }

    //get all the items in the screen
    QList<QGraphicsItem *> allItems = this->items();

    //remove all items from the scene and delete pointers
    for (QGraphicsItem *item : allItems) {
        // Remove the item from the scene
        this->removeItem(item);
        // Delete the item
        delete item;
        //qDebug() << "removed item";
    }


    //close the window
    if (this->mainWidget) {
        //qDebug() << "closing main widget";
        mainWidget->close();
        //qDebug() << "successfully closed";
    }
}

