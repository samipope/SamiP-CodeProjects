#include "controller.h"

//Constructor
Controller::Controller() {
    //default level
    level = 1;

    //initialize death screen
    gameOverScreen = new deathScreen();

    //Connect the signal from the login to a slot in the Controller widget
    connect(&loginScreen, &login::loginSignal, this, &Controller::handleLogin, Qt::QueuedConnection);

    //connect the replay button signal to the slot in the Controller
    connect(gameOverScreen, &deathScreen::replayClicked, this, &Controller::handleReplay, Qt::QueuedConnection);
}


// Brings up the login screen
void Controller::setUpLogin(){
    loginScreen.show();
}


// Slot to handle the login event
void Controller::handleLogin(const QString &data, int levelNum) {
    //initialize the player
    player = new User();
    level = levelNum;

    //if not a guest, load the user data
    if(data != "Guest"){
        //qDebug() << "Not a guest";

        //load player data from file
        player->fromJson(data);
        //hides login screen
        loginScreen.hide();
        //sends the login success message
        player->signInMessage();
        //hides sign in message (for some reason we needed to do it twice)
        loginScreen.hide();
    }


    //make a new game
    game = new GameScene();
    //connect the "game over signal" to the controller slot
    connect(game, &GameScene::gameOverSignal, this, &Controller::handleGameOver, Qt::QueuedConnection);

    //start the game
    game->startGameScene(*player, level);
    //hide the login screen
    loginScreen.hide();
}

//slot to handle the gameover event
void Controller::handleGameOver(int score){
    //qDebug() << "Game over inside controller";

    //store the user data before deleting the game
    User finalUser = game->gameUser;
    //takes care of the game window
    game->takeDownGame();

    //get rid of the game
    game->clear();
    delete game;
    game = nullptr;
    //qDebug() << "deleted game";

    //pull up the death screen
    gameOverScreen->displayDeath(score, finalUser);

    //qDebug() << "reached end of handle game over";
}

//slot to handle the replay event
void Controller::handleReplay(){
    //hide the gameover screen
    gameOverScreen->hide();

    //make a new game
    game = new GameScene();
    //connect the game over signal to the controller
    connect(game, &GameScene::gameOverSignal, this, &Controller::handleGameOver, Qt::QueuedConnection);

    //start the new game
    game->startGameScene(*player, level);
}
