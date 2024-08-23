#include "user.h"
#include <QMessageBox>

/**
 * @brief User::User default constructor makes guest?
 */
User::User() {
    firstName = "None";
    lastName = "None";
    userName = "Guest";
    password = "Nopassword1";
    birthDate = QDate(2000, 1, 1);
    gender = "female";
    profilePicPath = "://images/sami-hacker.png";
    QList<int> intList({0});
    scoreHistory = intList;
}

/**
 * @brief User::userAlreadyExists checks json hashmap to see if user already exists
 * @param possibleUserName string key
 * @return true if user already exists
 */
bool User::userAlreadyExists(QString possibleUserName) {
    // get dictionary from file
    QJsonObject loginObj = getJsonObjectFromFile(getPlayersDir() + "userslist.json");
    QMap<QString, QString> userMap = deserializeMap(loginObj);

    // check if key exists
    return userMap.contains(possibleUserName);
}

/**
 * @brief User::saveNewUser saves a new user's info in the userslist.json
 * @param newUserName
 * @param newPassword
 */
void User::saveNewUser(QString newUserName, QString newPassword) {
    // get dictionary from file
    QJsonObject loginObj = getJsonObjectFromFile(getPlayersDir() + "userslist.json");
    QMap<QString, QString> userMap = deserializeMap(loginObj);

    // add this user to the dict
    if(!userMap.contains(newUserName)) {
        userMap.insert(newUserName, newPassword);
    }

    // re-save to file
    QJsonObject newloginObj = serializeMap(userMap);
    saveJsonFile(getPlayersDir() + "userslist.json", newloginObj);

}

/**
 * @brief User::checkPassword checks if the provided password matches what's in the userslist.json. returns false
 * if the user is not found in the json
 * @param testUserName
 * @param testPassword
 * @return
 */
bool User::checkPassword(QString testUserName, QString testPassword) {
    // get dictionary from file
    QJsonObject loginObj = getJsonObjectFromFile(getPlayersDir() + "userslist.json");
    QMap<QString, QString> userMap = deserializeMap(loginObj);

    // add this user to the dict
    if(userMap.contains(testUserName)) {
        return testPassword == userMap[testUserName];
    }

    return false;
}


/**
 * @brief User::fromJson "constructor" to reset this User object using values from a json
 * @param playerUserName will find the json file for this user
 */
void User::fromJson(QString playerUserName) {

    QJsonObject userJson = getJsonObjectFromFile(getPlayersDir() + getHash(playerUserName) + ".json");

    // set identifiers
    firstName = userJson.value("firstName").toString();
    lastName = userJson.value("lastName").toString();
    userName = userJson.value("userName").toString();
    password = userJson.value("password").toString();

    // set birthdate
    birthDate = QDate::fromString(userJson.value("birthDate").toString(),Qt::TextDate);

    // set gender and profile pic path
    gender = userJson.value("gender").toString();
    profilePicPath = userJson.value("profilePicPath").toString();

    // fill in history
    QJsonArray scoreArray = userJson.value("scoreHistory").toArray();
    // Iterate over the array and add the integers to the QList
    scoreHistory.clear();
    for (const QJsonValue &value : scoreArray) {
        if (value.isDouble()) {
            scoreHistory.append(value.toInt());
        }
    }
}


/**
 * @brief User::UpdateScores compares the user's current score with the High Score  and updates high score if necessary.
 * Also updates the user's scores.
 */
void User::updateScores(int newScore){
    scoreHistory.append(newScore);

    // get current high score from file
    QString filePath = getPlayersDir() + "highscore.json";
    QJsonObject highScoreJson = getJsonObjectFromFile(filePath);
    int currentHighScore = highScoreJson.value("score").toInt();

    // if this is a new high score, update the file
    if(newScore > currentHighScore) {
        QJsonObject newHighScore;
        newHighScore["username"] = userName;
        newHighScore["score"] = newScore;

        saveJsonFile(filePath, newHighScore);
    }
}

/**
 * @brief User::toJson makes a QJsonObject representing this User
 * @return
 */
QJsonObject User::toJson()  {
    QJsonObject obj;

    // get basic fields
    obj["firstName"] = firstName;
    obj["lastName"] = lastName;
    obj["userName"] = userName;
    obj["password"] = password;
    obj["birthDate"] = birthDate.toString(Qt::TextDate);
    obj["gender"] = gender;
    obj["profilePicPath"] = profilePicPath;

    //jsonify the score list
    QJsonArray jsonArray;
    for (const int& score : scoreHistory) {
        jsonArray.append(score);
    }

    obj["scoreHistory"] = jsonArray;

    return obj;
}

/**
 * @brief User::save This save the User to a json file in the "Players" directory. The filename is a hex hash of the username
 */
void User::save() {
    if(userName == "Guest") {
        return;
    }
    // Save the JSON document to a file in the "Players/" directory
    QString filePath = getPlayersDir() + getHash(userName) + ".json";

    saveJsonFile(filePath, this->toJson());
}

/**
 * @brief User::saveJsonFile helper method to save a json object to a file
 * @param filePath path for the file
 * @param jsonObject
 */
void User::saveJsonFile(QString filePath, QJsonObject jsonObject) {
    QJsonDocument jsonDoc(jsonObject);

    QFile file(filePath);
    if (file.open(QIODevice::WriteOnly)) {
        file.write(jsonDoc.toJson());
        file.close();
    } else {
        qDebug() << "Failed to save JSON document";
    }
}

/**
 * @brief User::getPlayersDir Helper function to find the file path of the Players/ directory on this machine
 * @return QString of complete file path
 */
QString User::getPlayersDir() {
    // Get the full path of the current source file
    QString currentFile = QFileInfo(__FILE__).absoluteFilePath();

    // Navigate up the directory structure until we find the project file
    QDir projectDir = QDir(currentFile);
    do {
        projectDir.cdUp();
    } while (!projectDir.exists("WaterDropletGame.pro"));

    QString directoryPath = projectDir.path();
    QString usersDir = directoryPath + "/Players/";

    // Create the "Users/" directory if it doesn't exist
    QDir dir(usersDir);
    if (!dir.exists()) {
        dir.mkpath(".");
    }

    return usersDir;
}

/**
 * @brief User::getJsonObjectFromFile Returns a QJsonOject from a file if it is a correctly formatted Json file, otherwise
 * returns a blank object
 * @param filePath complete file path to json file
 * @return
 */
QJsonObject User::getJsonObjectFromFile(QString filePath) {

    QFile file(filePath);
    QJsonObject jsonObject;

    if (file.open(QIODevice::ReadOnly))
    {
        QByteArray content = file.readAll();
        file.close();

        QJsonDocument doc = QJsonDocument::fromJson(content); // use 2nd argument here to get the parsing error in case the input JSON is malformed
        jsonObject = doc.object();
    }  else {
        qDebug() << "couldn't open json file: " << filePath;
    }

    return jsonObject;

}

/**
 * @brief User::getHash A QString of a hex number representing the hash. This is created from the usernames,
 * assuming all Users have unique usernames
 * @param userNameString
 * @return
 */
QString User::getHash(QString userNameString) {
    QByteArray utf8Data = userNameString.toUtf8();

    // Compute the hash
    QByteArray hash = QCryptographicHash::hash(utf8Data, QCryptographicHash::Md5);

    // Convert the hash to a hexadecimal string
    QString hashString = QString(hash.toHex());
    return hashString;
}

/**
 * @brief User::serializeMap helper function to turn a hashmap to a json
 * @param dictionary
 * @return
 */
QJsonObject User::serializeMap(const QMap<QString, QString> &dictionary)
{
    QJsonObject obj;
    for (auto it = dictionary.constBegin(); it != dictionary.constEnd(); ++it) {
        obj.insert(it.key(), it.value());
    }
    return obj;
}

/**
 * @brief User::deserializeMap helper function to turn a json to a hashmap
 * @param obj
 * @return
 */
QMap<QString, QString> User::deserializeMap(const QJsonObject &obj)
{
    QMap<QString, QString> dictionary;
    for (auto it = obj.constBegin(); it != obj.constEnd(); ++it) {
        if (it.value().isString()) {
            dictionary.insert(it.key(), it.value().toString());
        }
    }
    return dictionary;
}




void User::signInMessage(){

    QDate currentDate = QDate::currentDate();

    if(currentDate.day() == birthDate.day() && currentDate.month() == birthDate.month()){
        QMessageBox birthdayMessage;

        //set window title
        birthdayMessage.setWindowTitle("Happy Birthday!");

        //set the message text
        birthdayMessage.setText("Happy Birthday " + userName);

        //set the icon
        birthdayMessage.setIconPixmap(QPixmap(":/images/birthday-cake.jpg"));

        //add buttons
        birthdayMessage.addButton(QMessageBox::Ok);

        //show the message box
        birthdayMessage.exec();

    }

    else{

        QMessageBox signInMessage;

        //set window title
        signInMessage.setWindowTitle("Login success");

        // Create a QString to store the output
        QString output;
        QTextStream stream(&output);

        // Append each element of the list to the string
        stream << "[";
        for (int i = 0; i < scoreHistory.size(); ++i) {
            stream << scoreHistory.at(i);
            if (i < scoreHistory.size() - 1) {
                stream << ", ";
            }
        }
        stream << "]";

        //set the message text
        signInMessage.setText("Welcome " + userName + "!!!\nYour Score History: " + output);

        //set the icon
        signInMessage.setIconPixmap(QPixmap(profilePicPath));

        //add buttons
        signInMessage.addButton(QMessageBox::Ok);

        //show the message box
        signInMessage.exec();

    }

}

//returns the user's high score
int User::getUserHighScore(){
    // Sorting the vector
    std::sort(scoreHistory.begin(), scoreHistory.end());

    return scoreHistory.last();
}

//returns the global high score
int User::getGlobalHighScore(){
    // get current high score from file
    QString filePath = getPlayersDir() + "highscore.json";
    QJsonObject highScoreJson = getJsonObjectFromFile(filePath);
    int currentHighScore = highScoreJson.value("score").toInt();

    return currentHighScore;
}
