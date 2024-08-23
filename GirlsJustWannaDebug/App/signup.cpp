#include "signup.h"
#include "ui_signup.h"
#include <QMessageBox>
#include <QFileDialog>
#include "user.h"
#include <QImage>
#include <QFile>
#include <QDir>

signup::signup(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::signup)
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



signup::~signup()
{
    delete ui;
}

bool isPasswordValid(const QString& password) {
    if (password.length() < 8) {
        return false;
    }

    bool hasLower = false, hasUpper = false, hasDigit = false;
    for (int i = 0; i < password.length(); ++i) {
        if (password[i].isLower()) hasLower = true;
        if (password[i].isUpper()) hasUpper = true;
        if (password[i].isDigit()) hasDigit = true;

        if (hasLower && hasUpper && hasDigit) {
            return true;
        }
    }

    return false;
}



void signup::on_OKButton_accepted()
{
    bool isUserValid = true;

    //TODO add a picture !!

    QString password = ui->passwordEntrySignup->text();
    if (!isPasswordValid(password)) {
        QMessageBox::warning(this, "Invalid Password", "Password must be at least 8 characters long and include at least one uppercase letter, one lowercase letter, and one digit.");
        return;
    }
    QString firstName = ui->firstNameEntry->text();
    if(firstName.isEmpty() || firstName.isNull()){
        QMessageBox::warning(this, "First name cannot be empty", "First name cannot be empty. Please try again");
        return;}


    QString lastName = ui->lastNameEntry->text();
    if(lastName.isEmpty() || lastName.isNull()){
        QMessageBox::warning(this, "Last name cannot be empty", "Last name cannot be empty. Please try again");
        return;
    }
    QDate dob = ui->DOBEntry->date();
    if(dob.isNull()){
        QMessageBox::warning(this, "DOB name cannot be empty", "Date of Birth cannot be empty. Please try again");
        return;
    }
     //how to do gender?

    QString gender = ui->genderBox->currentText();
    if(gender.isEmpty()) {
        QMessageBox::warning(this, "Gender not selected", "Gender not selected. Please select a gender");
        return;
    }

    QString userName = ui->usernameEntrySignup->text();
    if(userName.isEmpty() || userName.isNull()){
        QMessageBox::warning(this, "Username cannot be empty", "Username cannot be empty. Please try again");
        return;
    }
    //checking username not guest or Guest or GUEST
    else if(QString::compare(userName, "guest", Qt::CaseInsensitive) == 0){ //username cannot be guest because it messes with guest constructor
        QMessageBox::warning(this, "Username is invalid", "Username cannot be Guest or guest");
        return;
    }



    //CALL JSON function to see if user is valid here

    isUserValid = !User::userAlreadyExists(userName);

    if(isUserValid){ //if user is valid, call constructor
        // User *user = new User();
        // User newUser = new User(username, password, firstName, lastName, DOB, gender)
        User newUser;
        newUser.firstName = firstName;
        newUser.lastName = lastName;
        newUser.userName = userName;
        newUser.password = password;
        newUser.birthDate = dob;
        newUser.gender = gender;

        //save profile pic to players Directory and write path
        QString fileName = User::getHash(newUser.userName);
        fileName += profilePicExtension;
        newUser.profilePicPath = saveImage(profilePic, fileName);

        //save the user
        newUser.save();
        newUser.saveNewUser(newUser.userName, newUser.password);
        QMessageBox::information(this, "Registration Complete", "Thanks for signing up!");

    }
    else{
        QMessageBox::warning(this, "Username already taken", "This username is taken, please choose another");
        return;
    }



    accept();

}


void signup::on_pushButton_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this,tr("choose"),"",tr("Images(*.png *.jpg *.jpeg * .gif)"));
    if(QString::compare(filename,QString())!=0){
        //QImage image;
        profilePicExtension = getFileExtension(filename);
        profilePic.load(filename);

    }
}

// saves an image to the player directory and returns the full path to the file
QString saveImage(const QImage &newimage, const QString &fileName){
    //construct the file path
    QString folderPath = User::getPlayersDir();
    QString filePath = folderPath + fileName;

    //Save the image to the file
    bool saved = newimage.save(filePath);

    if(saved){
        qDebug() << "Image saved successfully to: " << filePath;
    }
    else{
        qDebug() << "Failed to save image to:" << filePath;
    }
    return filePath;
}


//returns the extension of the given file name
QString getFileExtension(const QString &fileName){
    QStringList list = fileName.split(".");

    return "." + list.last();
}

