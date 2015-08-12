/********************************************************************************
** Form generated from reading UI file 'gui.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GUI_H
#define UI_GUI_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QPushButton>
#include <QtGui/QStatusBar>
#include <QtGui/QToolButton>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Gui
{
public:
    QWidget *centralwidget;
    QGridLayout *gridLayout;
    QLabel *label;
    QLineEdit *isoEd;
    QLabel *label_2;
    QLineEdit *fragsEd;
    QLabel *label_3;
    QLineEdit *thetaEd;
    QLabel *label_4;
    QLineEdit *rhoEd;
    QLabel *label_5;
    QLineEdit *deltaEd;
    QLabel *label_10;
    QLineEdit *popSizeEd;
    QToolButton *popSizeShow;
    QFrame *line;
    QLabel *label_6;
    QLineEdit *dataTo;
    QToolButton *dataSelect;
    QLabel *label_9;
    QLineEdit *ltTo;
    QToolButton *ltSelect;
    QLabel *label_7;
    QLineEdit *cgTo;
    QToolButton *cgSelect;
    QLabel *label_8;
    QLineEdit *dotTo;
    QToolButton *dotSelect;
    QCheckBox *amEd;
    QHBoxLayout *hboxLayout;
    QPushButton *help;
    QPushButton *run;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *Gui)
    {
        if (Gui->objectName().isEmpty())
            Gui->setObjectName(QString::fromUtf8("Gui"));
        Gui->resize(389, 443);
        centralwidget = new QWidget(Gui);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout = new QGridLayout(centralwidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label = new QLabel(centralwidget);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        isoEd = new QLineEdit(centralwidget);
        isoEd->setObjectName(QString::fromUtf8("isoEd"));
        isoEd->setFrame(true);
        isoEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(isoEd, 0, 1, 1, 3);

        label_2 = new QLabel(centralwidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        fragsEd = new QLineEdit(centralwidget);
        fragsEd->setObjectName(QString::fromUtf8("fragsEd"));
        fragsEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(fragsEd, 1, 1, 1, 3);

        label_3 = new QLabel(centralwidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout->addWidget(label_3, 2, 0, 1, 1);

        thetaEd = new QLineEdit(centralwidget);
        thetaEd->setObjectName(QString::fromUtf8("thetaEd"));
        thetaEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(thetaEd, 2, 1, 1, 3);

        label_4 = new QLabel(centralwidget);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout->addWidget(label_4, 3, 0, 1, 1);

        rhoEd = new QLineEdit(centralwidget);
        rhoEd->setObjectName(QString::fromUtf8("rhoEd"));
        rhoEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(rhoEd, 3, 1, 1, 3);

        label_5 = new QLabel(centralwidget);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout->addWidget(label_5, 4, 0, 1, 1);

        deltaEd = new QLineEdit(centralwidget);
        deltaEd->setObjectName(QString::fromUtf8("deltaEd"));
        deltaEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(deltaEd, 4, 1, 1, 3);

        label_10 = new QLabel(centralwidget);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        gridLayout->addWidget(label_10, 5, 0, 1, 1);

        popSizeEd = new QLineEdit(centralwidget);
        popSizeEd->setObjectName(QString::fromUtf8("popSizeEd"));
        popSizeEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(popSizeEd, 5, 1, 1, 1);

        popSizeShow = new QToolButton(centralwidget);
        popSizeShow->setObjectName(QString::fromUtf8("popSizeShow"));

        gridLayout->addWidget(popSizeShow, 5, 2, 1, 2);

        line = new QFrame(centralwidget);
        line->setObjectName(QString::fromUtf8("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        gridLayout->addWidget(line, 6, 0, 1, 4);

        label_6 = new QLabel(centralwidget);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout->addWidget(label_6, 7, 0, 1, 1);

        dataTo = new QLineEdit(centralwidget);
        dataTo->setObjectName(QString::fromUtf8("dataTo"));

        gridLayout->addWidget(dataTo, 7, 1, 1, 2);

        dataSelect = new QToolButton(centralwidget);
        dataSelect->setObjectName(QString::fromUtf8("dataSelect"));

        gridLayout->addWidget(dataSelect, 7, 3, 1, 1);

        label_9 = new QLabel(centralwidget);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        gridLayout->addWidget(label_9, 8, 0, 1, 1);

        ltTo = new QLineEdit(centralwidget);
        ltTo->setObjectName(QString::fromUtf8("ltTo"));

        gridLayout->addWidget(ltTo, 8, 1, 1, 2);

        ltSelect = new QToolButton(centralwidget);
        ltSelect->setObjectName(QString::fromUtf8("ltSelect"));

        gridLayout->addWidget(ltSelect, 8, 3, 1, 1);

        label_7 = new QLabel(centralwidget);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        gridLayout->addWidget(label_7, 9, 0, 1, 1);

        cgTo = new QLineEdit(centralwidget);
        cgTo->setObjectName(QString::fromUtf8("cgTo"));

        gridLayout->addWidget(cgTo, 9, 1, 1, 2);

        cgSelect = new QToolButton(centralwidget);
        cgSelect->setObjectName(QString::fromUtf8("cgSelect"));

        gridLayout->addWidget(cgSelect, 9, 3, 1, 1);

        label_8 = new QLabel(centralwidget);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout->addWidget(label_8, 10, 0, 1, 1);

        dotTo = new QLineEdit(centralwidget);
        dotTo->setObjectName(QString::fromUtf8("dotTo"));

        gridLayout->addWidget(dotTo, 10, 1, 1, 2);

        dotSelect = new QToolButton(centralwidget);
        dotSelect->setObjectName(QString::fromUtf8("dotSelect"));

        gridLayout->addWidget(dotSelect, 10, 3, 1, 1);

        amEd = new QCheckBox(centralwidget);
        amEd->setObjectName(QString::fromUtf8("amEd"));

        gridLayout->addWidget(amEd, 11, 0, 1, 3);

        hboxLayout = new QHBoxLayout();
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        help = new QPushButton(centralwidget);
        help->setObjectName(QString::fromUtf8("help"));

        hboxLayout->addWidget(help);

        run = new QPushButton(centralwidget);
        run->setObjectName(QString::fromUtf8("run"));

        hboxLayout->addWidget(run);


        gridLayout->addLayout(hboxLayout, 12, 0, 1, 4);

        Gui->setCentralWidget(centralwidget);
        statusbar = new QStatusBar(Gui);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        Gui->setStatusBar(statusbar);

        retranslateUi(Gui);

        QMetaObject::connectSlotsByName(Gui);
    } // setupUi

    void retranslateUi(QMainWindow *Gui)
    {
        Gui->setWindowTitle(QApplication::translate("Gui", "SimMLST", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("Gui", "Number of isolates", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_WHATSTHIS
        isoEd->setWhatsThis(QString());
#endif // QT_NO_WHATSTHIS
        isoEd->setInputMask(QString());
        isoEd->setText(QApplication::translate("Gui", "100", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("Gui", "Lengths of gene fragments", 0, QApplication::UnicodeUTF8));
        fragsEd->setText(QApplication::translate("Gui", "400,400,400,400,400,400,400", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("Gui", "Mutation rate (theta)", 0, QApplication::UnicodeUTF8));
        thetaEd->setText(QApplication::translate("Gui", "100", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("Gui", "Recombination rate (rho)", 0, QApplication::UnicodeUTF8));
        rhoEd->setText(QApplication::translate("Gui", "100", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("Gui", "Mean tract length (delta)", 0, QApplication::UnicodeUTF8));
        deltaEd->setText(QApplication::translate("Gui", "500", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("Gui", "Population size model", 0, QApplication::UnicodeUTF8));
        popSizeShow->setText(QApplication::translate("Gui", "Show", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("Gui", "Save data to...", 0, QApplication::UnicodeUTF8));
        dataSelect->setText(QApplication::translate("Gui", "...", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("Gui", "Save local trees to...", 0, QApplication::UnicodeUTF8));
        ltSelect->setText(QApplication::translate("Gui", "...", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("Gui", "Save clonal genealogy to...", 0, QApplication::UnicodeUTF8));
        cgSelect->setText(QApplication::translate("Gui", "...", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("Gui", "Save DOT graph to...", 0, QApplication::UnicodeUTF8));
        dotSelect->setText(QApplication::translate("Gui", "...", 0, QApplication::UnicodeUTF8));
        amEd->setText(QApplication::translate("Gui", "Save ancestral material in the DOT graph", 0, QApplication::UnicodeUTF8));
        help->setText(QApplication::translate("Gui", "Help", 0, QApplication::UnicodeUTF8));
        run->setText(QApplication::translate("Gui", "Simulate!", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Gui: public Ui_Gui {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GUI_H
