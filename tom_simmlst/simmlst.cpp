#include <QApplication>
#include "guiimpl.h"
#include "arg.h"
#include "rng.h"

static const char * help=
    "\
    Usage: simMLST [OPTIONS]\n\
    \n\
    Options:\n\
    -N NUM   Sets the number of isolates (default is 100)\n\
    -T NUM   Sets the value of theta (default is 100)\n\
    -m NUM   Sets the minimum probability of mutation in an interval of external recombination between 0 & 1 (default is 0)\n\
    -M NUM   Sets the maximum probability of mutation in an interval of external recombination between 0 & 1 (default is 0)\n\
    -R NUM   Sets the value of rho (default is 100)\n\
    -r NUM   Sets the rate of external recombination (default is 0)\n\
    -D NUM   Sets the value of delta (default is 500)\n\
    -e NUM   Sets the average length of external recombinant interval (default is 0)\n\
    -B NUM,...,NUM Sets the number and length of the fragments\n\
             (default is 10000)\n\
    -G NUM   Sets the gap between each fragment(default is 0)\n\
    -C T,N   Sets the population size constant and equal to N \n\
             further in the past than time T.  Negative values of N\n\
             are interpreted as the current population size\n\
    -E T,R   Sets the population size exponentially growing with\n\
             rate R further in the past than time T\n\
    -s NUM   Use given seed to initiate random number generator\n\
    -o FILE  Export data to given file\n\
    -c FILE  Export clonal genealogy to given file\n\
    -l FILE  Export local trees to given file\n\
    -b FILE  Write log file of recombinant break interval locations\n\
    -d FILE  Export DOT graph to given file\n\
    -a       Include ancestral material in the DOT graph\n\
    ";

int main(int argc, char *argv[]) {
    if (argc==1) {//GUI mode
    makerng();
    QApplication app(argc, argv);
    GuiImpl guiimpl;
    guiimpl.show();
    return app.exec();
    }

    //Initialization
    PopSize * popsize=NULL;
    int n=100;
    double theta=100.0;
    double theta_extMin=0.5;
    double theta_extMax=1.0;
    double rho=100.0;
    double rho_ext=0.0;
    double delta=500.0;
    double delta_ext=500.0;
    int seed=-1;
    bool am=false;
    vector<int> blocks;
    for (int i=0;i<2;i++) blocks.push_back(i*10000);
    vector<int> gaps;
    int c;
    string cfile="";//File to which the clonal genealogy is exported
    string lfile="";//File to which the local trees are exported
    string dfile="";//File to which the DOT graph is exported
    string ofile="";//File to which the data is exported
    string bfile="";//File to which recombinant intervals are exported
    double p1,p2;
    char * pch;
    while ((c = getopt (argc, argv, "ahN:T:m:M:R:r:D:e:s:B:G:c:l:d:o:C:E:")) != -1)
    switch (c)
    {
        case('N'):n=atoi(optarg);break;
        case('m'):theta_extMin=atof(optarg);break;
        case('M'):theta_extMax=atof(optarg);break;
        case('T'):theta=atof(optarg);break;
        case('R'):rho=atof(optarg);break;
        case('r'):rho_ext=atof(optarg);break;
        case('D'):delta=atof(optarg);break;
        case('e'):delta_ext=atof(optarg);break;
        case('s'):seed=atoi(optarg);break;
        case('B'):blocks=Arg::makeBlocks(optarg);break;
        case('G'):gaps=Arg::makeGaps(optarg);break;
        case('h'):cout<<help<<endl;return 1;break;
        case('a'):am=true;break;
        case('c'):cfile=optarg;break;
        case('l'):lfile=optarg;break;
        case('d'):dfile=optarg;break;
        case('o'):ofile=optarg;break;
        case('b'):bfile=optarg;break;
        case('C'):if (popsize==NULL) popsize=new PopSize();p1=atof(optarg);pch=strtok(optarg,",");pch=strtok(NULL,",");p2=atof(pch);popsize->addEvent(false,p1,p2);break;
        case('E'):if (popsize==NULL) popsize=new PopSize();p1=atof(optarg);pch=strtok(optarg,",");pch=strtok(NULL,",");p2=atof(pch);popsize->addEvent(true ,p1,p2);break;
        case '?':cout<<"Wrong arguments: did not recognise "<<c<<" "<<optarg<<endl<<help<<endl;return 1;
        default:abort();
    }
    if (seed==-1) makerng(); else {rng=gsl_rng_alloc(gsl_rng_default);gsl_rng_set(rng,seed);};
    if (gaps.size() == 0){
      for (size_t i=0;i<blocks.size()-1;++i) gaps.push_back(0);
    }

    if (gaps.size() != blocks.size()-1) {cout << "Wrong number of gaps given (must be same as number of blocks)" << endl;return 1;}

    if (argc-optind>0) {cout<<"Wrong arguments."<<endl<<help<<endl;return 1;}
    //Build the ARG
    Arg * arg=new Arg(n,rho,rho_ext,delta,delta_ext,blocks,gaps,popsize);
    //Build the data and export it
    if (ofile.length()>0) {
    Data * data=arg->drawData(theta,theta_extMin,theta_extMax);
    ofstream dat;
    dat.open(ofile.data());
    data->output(&dat);
    dat.close();
    delete(data);}
    //Extract the clonal genealogy and export it
    if (cfile.length()>0) {
    string truth=arg->extractCG();
    ofstream tru;
    tru.open(cfile.data());
    tru<<truth<<endl;
    tru.close();}
    //Extract the local trees and export them
    if (lfile.length()>0) {
    ofstream lf;
    lf.open(lfile.data());
    arg->outputLOCAL(&lf);
    lf.close();}
    //Export to DOT format
    if (dfile.length()>0) {
    ofstream dot;
    dot.open(dfile.data());
    arg->outputDOT(&dot,am);
    dot.close();}
    //Write recombinant breaks to file
    if (bfile.length()>0) {
    ofstream breaks;
    breaks.open(bfile.data());
    arg->outputBREAKS(&breaks);
    breaks.close();}
    delete(arg);
    if (popsize!=NULL) delete(popsize);
    gsl_rng_free(rng);
}
