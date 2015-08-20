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
    QLabel *label_20;
    QLineEdit *gapsEd;
    QLabel *label_3;
    QLineEdit *thetaEd;
    QLabel *label_31;
    QLineEdit *theta_extMinEd;
    QLabel *label_32;
    QLineEdit *theta_extMaxEd;
    QLabel *label_4;
    QLineEdit *rhoEd;
    QLabel *label_41;
    QLineEdit *rho_extEd;
    QLabel *label_5;
    QLineEdit *deltaEd;
    QLabel *label_51;
    QLineEdit *delta_extEd;
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
    QLabel *label_71;
    QLineEdit *rbTo;
    QToolButton *rbSelect;
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
        Gui->resize(443, 543);
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

        label_20 = new QLabel(centralwidget);
        label_20->setObjectName(QString::fromUtf8("label_20"));

        gridLayout->addWidget(label_20, 2, 0, 1, 1);

        gapsEd = new QLineEdit(centralwidget);
        gapsEd->setObjectName(QString::fromUtf8("gapsEd"));
        gapsEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(gapsEd, 2, 1, 1, 3);

        label_3 = new QLabel(centralwidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout->addWidget(label_3, 3, 0, 1, 1);

        thetaEd = new QLineEdit(centralwidget);
        thetaEd->setObjectName(QString::fromUtf8("thetaEd"));
        thetaEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(thetaEd, 3, 1, 1, 3);

        label_31 = new QLabel(centralwidget);
        label_31->setObjectName(QString::fromUtf8("label_31"));

        gridLayout->addWidget(label_31, 4, 0, 1, 1);

        theta_extMinEd = new QLineEdit(centralwidget);
        theta_extMinEd->setObjectName(QString::fromUtf8("theta_extMinEd"));
        theta_extMinEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(theta_extMinEd, 4, 1, 1, 3);

        label_32 = new QLabel(centralwidget);
        label_32->setObjectName(QString::fromUtf8("label_32"));

        gridLayout->addWidget(label_32, 5, 0, 1, 1);

        theta_extMaxEd = new QLineEdit(centralwidget);
        theta_extMaxEd->setObjectName(QString::fromUtf8("theta_extMaxEd"));
        theta_extMaxEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(theta_extMaxEd, 5, 1, 1, 3);

        label_4 = new QLabel(centralwidget);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout->addWidget(label_4, 6, 0, 1, 1);

        rhoEd = new QLineEdit(centralwidget);
        rhoEd->setObjectName(QString::fromUtf8("rhoEd"));
        rhoEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(rhoEd, 6, 1, 1, 3);

        label_41 = new QLabel(centralwidget);
        label_41->setObjectName(QString::fromUtf8("label_41"));

        gridLayout->addWidget(label_41, 7, 0, 1, 1);

        rho_extEd = new QLineEdit(centralwidget);
        rho_extEd->setObjectName(QString::fromUtf8("rho_extEd"));
        rho_extEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(rho_extEd, 7, 1, 1, 3);

        label_5 = new QLabel(centralwidget);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout->addWidget(label_5, 8, 0, 1, 1);

        deltaEd = new QLineEdit(centralwidget);
        deltaEd->setObjectName(QString::fromUtf8("deltaEd"));
        deltaEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(deltaEd, 8, 1, 1, 3);

        label_51 = new QLabel(centralwidget);
        label_51->setObjectName(QString::fromUtf8("label_51"));

        gridLayout->addWidget(label_51, 9, 0, 1, 1);

        delta_extEd = new QLineEdit(centralwidget);
        delta_extEd->setObjectName(QString::fromUtf8("delta_extEd"));
        delta_extEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(delta_extEd, 9, 1, 1, 3);

        label_10 = new QLabel(centralwidget);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        gridLayout->addWidget(label_10, 10, 0, 1, 1);

        popSizeEd = new QLineEdit(centralwidget);
        popSizeEd->setObjectName(QString::fromUtf8("popSizeEd"));
        popSizeEd->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(popSizeEd, 10, 1, 1, 1);

        popSizeShow = new QToolButton(centralwidget);
        popSizeShow->setObjectName(QString::fromUtf8("popSizeShow"));

        gridLayout->addWidget(popSizeShow, 10, 2, 1, 2);

        line = new QFrame(centralwidget);
        line->setObjectName(QString::fromUtf8("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        gridLayout->addWidget(line, 11, 0, 1, 4);

        label_6 = new QLabel(centralwidget);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout->addWidget(label_6, 12, 0, 1, 1);

        dataTo = new QLineEdit(centralwidget);
        dataTo->setObjectName(QString::fromUtf8("dataTo"));

        gridLayout->addWidget(dataTo, 12, 1, 1, 2);

        dataSelect = new QToolButton(centralwidget);
        dataSelect->setObjectName(QString::fromUtf8("dataSelect"));

        gridLayout->addWidget(dataSelect, 12, 3, 1, 1);

        label_9 = new QLabel(centralwidget);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        gridLayout->addWidget(label_9, 13, 0, 1, 1);

        ltTo = new QLineEdit(centralwidget);
        ltTo->setObjectName(QString::fromUtf8("ltTo"));

        gridLayout->addWidget(ltTo, 13, 1, 1, 2);

        ltSelect = new QToolButton(centralwidget);
        ltSelect->setObjectName(QString::fromUtf8("ltSelect"));

        gridLayout->addWidget(ltSelect, 13, 3, 1, 1);

        label_7 = new QLabel(centralwidget);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        gridLayout->addWidget(label_7, 14, 0, 1, 1);

        cgTo = new QLineEdit(centralwidget);
        cgTo->setObjectName(QString::fromUtf8("cgTo"));

        gridLayout->addWidget(cgTo, 14, 1, 1, 2);

        cgSelect = new QToolButton(centralwidget);
        cgSelect->setObjectName(QString::fromUtf8("cgSelect"));

        gridLayout->addWidget(cgSelect, 14, 3, 1, 1);

        label_71 = new QLabel(centralwidget);
        label_71->setObjectName(QString::fromUtf8("label_71"));

        gridLayout->addWidget(label_71, 15, 0, 1, 1);

        rbTo = new QLineEdit(centralwidget);
        rbTo->setObjectName(QString::fromUtf8("rbTo"));

        gridLayout->addWidget(rbTo, 15, 1, 1, 2);

        rbSelect = new QToolButton(centralwidget);
        rbSelect->setObjectName(QString::fromUtf8("rbSelect"));

        gridLayout->addWidget(rbSelect, 15, 3, 1, 1);

        label_8 = new QLabel(centralwidget);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout->addWidget(label_8, 16, 0, 1, 1);

        dotTo = new QLineEdit(centralwidget);
        dotTo->setObjectName(QString::fromUtf8("dotTo"));

        gridLayout->addWidget(dotTo, 16, 1, 1, 2);

        dotSelect = new QToolButton(centralwidget);
        dotSelect->setObjectName(QString::fromUtf8("dotSelect"));

        gridLayout->addWidget(dotSelect, 16, 3, 1, 1);

        amEd = new QCheckBox(centralwidget);
        amEd->setObjectName(QString::fromUtf8("amEd"));

        gridLayout->addWidget(amEd, 17, 0, 1, 3);

        hboxLayout = new QHBoxLayout();
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        help = new QPushButton(centralwidget);
        help->setObjectName(QString::fromUtf8("help"));

        hboxLayout->addWidget(help);

        run = new QPushButton(centralwidget);
        run->setObjectName(QString::fromUtf8("run"));

        hboxLayout->addWidget(run);


        gridLayout->addLayout(hboxLayout, 18, 0, 1, 4);

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
        fragsEd->setText(QApplication::translate("Gui", "10000", 0, QApplication::UnicodeUTF8));
        label_20->setText(QApplication::translate("Gui", "Gap sizes (must be same number as fragments)", 0, QApplication::UnicodeUTF8));
        gapsEd->setText(QApplication::translate("Gui", "0", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("Gui", "Mutation rate (theta)", 0, QApplication::UnicodeUTF8));
        thetaEd->setText(QApplication::translate("Gui", "100", 0, QApplication::UnicodeUTF8));
        label_31->setText(QApplication::translate("Gui", "External minimum mutation probability (theta_extMin)", 0, QApplication::UnicodeUTF8));
        theta_extMinEd->setText(QApplication::translate("Gui", "0.5", 0, QApplication::UnicodeUTF8));
        label_32->setText(QApplication::translate("Gui", "External maximum mutation probability (theta_extMax)", 0, QApplication::UnicodeUTF8));
        theta_extMaxEd->setText(QApplication::translate("Gui", "1", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("Gui", "Recombination rate (rho)", 0, QApplication::UnicodeUTF8));
        rhoEd->setText(QApplication::translate("Gui", "100", 0, QApplication::UnicodeUTF8));
        label_41->setText(QApplication::translate("Gui", "External recombination rate (rho_ext)", 0, QApplication::UnicodeUTF8));
        rho_extEd->setText(QApplication::translate("Gui", "0", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("Gui", "Mean tract length (delta)", 0, QApplication::UnicodeUTF8));
        deltaEd->setText(QApplication::translate("Gui", "500", 0, QApplication::UnicodeUTF8));
        label_51->setText(QApplication::translate("Gui", "External mean tract length (delta_ext)", 0, QApplication::UnicodeUTF8));
        delta_extEd->setText(QApplication::translate("Gui", "0", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("Gui", "Population size model", 0, QApplication::UnicodeUTF8));
        popSizeShow->setText(QApplication::translate("Gui", "Show", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("Gui", "Save data to...", 0, QApplication::UnicodeUTF8));
        dataSelect->setText(QApplication::translate("Gui", "...", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("Gui", "Save local trees to...", 0, QApplication::UnicodeUTF8));
        ltSelect->setText(QApplication::translate("Gui", "...", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("Gui", "Save clonal genealogy to...", 0, QApplication::UnicodeUTF8));
        cgSelect->setText(QApplication::translate("Gui", "...", 0, QApplication::UnicodeUTF8));
        label_71->setText(QApplication::translate("Gui", "Save recombinant intervals to...", 0, QApplication::UnicodeUTF8));
        rbSelect->setText(QApplication::translate("Gui", "...", 0, QApplication::UnicodeUTF8));
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
