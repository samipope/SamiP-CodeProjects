#ifndef SIGNUP_H
#define SIGNUP_H

#include <QDialog>
#include <QImage>

namespace Ui {
class signup;
}

class signup : public QDialog
{
    Q_OBJECT

public:
    explicit signup(QWidget *parent = nullptr);
    ~signup();

private:
    QImage profilePic;
    QString profilePicExtension;

private slots:

    void on_OKButton_accepted();

    void on_pushButton_clicked();

private:
    Ui::signup *ui;
};

QString saveImage(const QImage &newimage, const QString &fileName);
QString getFileExtension(const QString &fileName);


#endif // SIGNUP_H
