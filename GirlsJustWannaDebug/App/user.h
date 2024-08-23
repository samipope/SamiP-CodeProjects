#ifndef USER_H
#define USER_H

#include <QCoreApplication>
#include <QImage>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QByteArray>
#include <QDebug>
#include <QJsonArray>
#include <QBuffer>
#include <QFile>
#include <QDir>
#include <QFileInfo>
#include <QCryptographicHash>
#include <QStandardPaths>

class User
{
private:
    static QJsonObject getJsonObjectFromFile(QString filePath);
    static QJsonObject serializeMap(const QMap<QString, QString> &dictionary);
    static QMap<QString, QString> deserializeMap(const QJsonObject &obj);
    static void saveJsonFile(QString filePath, QJsonObject jsonObject);
    QJsonObject toJson();

public:
    QString firstName;
    QString lastName;
    QString userName;
    QString password;
    QDate birthDate;
    QString gender;
    QString profilePicPath;
    QList<int> scoreHistory;


    User();
    void fromJson(QString playerUserName);
    void updateScores(int newScore);
    static User* getGuest();
    void save();
    static bool userAlreadyExists(QString possibleUserName);
    static bool checkPassword(QString testUserName, QString testPassword);
    static void saveNewUser(QString newUserName, QString newPassword);
    static QString getPlayersDir();
    static QString getHash(QString userNameString);
    void signInMessage();
    int getUserHighScore();
    static int getGlobalHighScore();
};

#endif // USER_H
