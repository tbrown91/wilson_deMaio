#include "guiimpl.h"
//
GuiImpl::GuiImpl( QWidget * parent, Qt::WFlags f)
    : QMainWindow(parent, f) {
  setupUi(this);
  isoEd->setValidator(new QIntValidator(2, 999999999, this));
  thetaEd->setValidator(new QDoubleValidator(0.0, 9999999999.0,10, this));
  theta_extMinEd->setValidator(new QDoubleValidator(0.0, 9999999999.0,10, this));
  theta_extMaxEd->setValidator(new QDoubleValidator(0.0, 9999999999.0,10, this));
  rhoEd->setValidator(  new QDoubleValidator(0.0, 9999999999.0,10, this));
  rho_extEd->setValidator(  new QDoubleValidator(0.0, 9999999999.0,10, this));
  deltaEd->setValidator(new QDoubleValidator(0.0, 9999999999.0,10, this));
  delta_extEd->setValidator(new QDoubleValidator(0.0, 9999999999.0,10, this));
  QRegExp rx("[0-9,]*");
  fragsEd->setValidator(new QRegExpValidator(rx,this));
  gapsEd->setValidator(new QRegExpValidator(rx,this));
}
//

void GuiImpl::on_help_clicked() {
  QMessageBox::about(0,"Help","\n\
Number of isolates:\tSets the number of isolates (default is 100)\n\
Length of gene fragments:\tSets the number and length of the fragments\n\
\t(default is 10000)\n\
Gap size:\tSets the gap size between gene fragements, must be same number as number of gene fragments:\t(default is 0)\n\
Mutation rate:\tSets the value of theta (default is 100)\n\
External minimum mutation probability:\tSets the minimum value of theta_ext (default is 0)\n\
External maximum mutation probability:\tSets the maximum value of theta_ext (default is 0)\n\
Recombination rate:\tSets the value of rho (default is 100)\n\
External recombination rate:\tSets the value of rho_ext (default is 0)\n\
Mean tract length delta:\tSets the value of delta (default is 500)\n\
Mean external tract length delta_ext:\tSets the value of delta_ext (default is 0)\n\
Population size model:\n\
\tThe population size parameters -C and -E are considered\n\
\tas events with only the most recent is active.  For example:\n\
\t-E 1,1 -C 2,-1\n\
\thas constant size from time 0 to 1, then exponential growth at\n\
\trate 1 from time 1 to 2, then constant size from time 2.\n\
\t-C T,N\tSets the population size constant and equal to N \n\
\tfurther in the past than time T.  Negative values of N\n\
\tare interpreted as the current population size\n\
\t-E T,R\tSets the population size exponentially growing with\n\
\trate R further in the past than time T\n\
\n\
Save Data to:\tExport sequence data to given file\n\
Save local trees to:\tExport local trees to given file\n\
Save clonal genealogy to:\tExport clonal genealogy to given file\n\
Save recombinatn intervals to:\tExport recombinant intervals to given file\n\
Save DOT graph to:\tExport DOT graph to given file\n\
Save ancestral material in the DOT graph:\n\
\tInclude ancestral material bars in the DOT graph.\n\
\tDisable this for large graphs to keep filesizes down\n\
\n\
");
}

void GuiImpl::on_dataSelect_clicked() {
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save data to file"),".","XMFA file (*.xmfa);;All files (*)");
  if (qstr==NULL) return;
  dataTo->setText(qstr);
}

void GuiImpl::on_cgSelect_clicked() {
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save genealogy to file"),".","Newick file (*.nwk);;All files (*)");
  if (qstr==NULL) return;
  cgTo->setText(qstr);
}

void GuiImpl::on_dotSelect_clicked() {
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save graph to file"),".","DOT file (*.dot);;All files (*)");
  if (qstr==NULL) return;
  dotTo->setText(qstr);
}

void GuiImpl::on_rbSelect_clicked() {
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save recombinant interval to file"),".","recomb file (*.log);;All files (*)");
  if (qstr==NULL) return;
  dotTo->setText(qstr);
}

void GuiImpl::on_ltSelect_clicked() {
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save local trees to file"),".","Newick file (*.nwk);;All files (*)");
  if (qstr==NULL) return;
  ltTo->setText(qstr);
}

void GuiImpl::on_run_clicked() {
  run->setText("Running...");run->setDisabled(true);run->repaint();QApplication::processEvents();
  int n=isoEd->text().toInt();
  double theta=thetaEd->text().toDouble();
  double theta_extMin=theta_extMinEd->text().toDouble();
  double theta_extMax=theta_extMaxEd->text().toDouble();
  double rho=rhoEd->text().toDouble();
  double rho_ext=rho_extEd->text().toDouble();
  double delta=deltaEd->text().toDouble();
  double delta_ext=delta_extEd->text().toDouble();
  vector<int> blocks=Arg::makeBlocks(fragsEd->text().toStdString());
  vector<int> gaps=Arg::makeBlocks(gapsEd->text().toStdString());
  PopSize * popsize=NULL;
  QString popstr=popSizeEd->text();
  //Interpret population size model
  if (!popstr.isEmpty()) {
      popsize=new PopSize();
      QRegExp qre("\\-([CE])\\ +([\\-0-9\\.]+),([\\-0-9\\.]+)");
      int cur=0;
      while (1) {cur=popstr.indexOf(qre,cur);if (cur<0) break;
          popsize->addEvent(qre.cap(1).compare("E")==0,qre.cap(2).toDouble(),qre.cap(3).toDouble());
          cur++;
        }
    }
  //if (popsize!=NULL) popsize->show();
  //Build the ARG
  Arg * arg=new Arg(n,rho,rho_ext,delta,delta_ext,blocks,gaps,popsize);
  //Build the data and export it
  if (!dataTo->text().isEmpty()) {
      Data * data=arg->drawData(theta,theta_extMin,theta_extMax);
      ofstream dat;
      dat.open(dataTo->text().toStdString().data());
      data->output(&dat);
      dat.close();
      delete(data);
    }
  //Extract the local trees and export them
  if (!ltTo->text().isEmpty()) {
      ofstream lf;
      lf.open(ltTo->text().toStdString().data());
      arg->outputLOCAL(&lf);
      lf.close();
    }
  //Extract the clonal genealogy and export it
  if (!cgTo->text().isEmpty()) {
      string truth=arg->extractCG();
      ofstream tru;
      tru.open(cgTo->text().toStdString().data());
      tru<<truth<<endl;
      tru.close();
    }
  //Extract the recombinant intervals and export them
  if (!rbTo->text().isEmpty()) {
      ofstream breaks;
      breaks.open(rbTo->text().toStdString().data());
      arg->outputBREAKS(&breaks);
      breaks.close();
    }
  //Export to DOT format
  if (!dotTo->text().isEmpty()) {
      ofstream dot;
      dot.open(dotTo->text().toStdString().data());
      arg->outputDOT(&dot,amEd->isChecked());
      dot.close();
    }
  delete(arg);
  run->setText("Simulate!");run->setDisabled(false);
  QMessageBox::about(0,"Information","All done!");
}

void GuiImpl::on_popSizeShow_clicked()
{
  //Interpret population size model
  QString popstr=popSizeEd->text();
  PopSize * ps=new PopSize();
  QRegExp qre("\\-([CE])\\ +([\\-0-9\\.]+),([\\-0-9\\.]+)");
  int cur=0;
  while (1) {cur=popstr.indexOf(qre,cur);if (cur<0) break;
      ps->addEvent(qre.cap(1).compare("E")==0,qre.cap(2).toDouble(),qre.cap(3).toDouble());
      cur++;
    }
  //Show population size model
  PlotImpl * pi=new PlotImpl();
  double maxtime;
  pi->setValues(ps->getVectorSizes(&maxtime));
  pi->setMaxX(maxtime);
  pi->show();
}
