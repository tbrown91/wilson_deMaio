#include "arg.h"
#include <ctime>
#include "coal_event.h"
#include "recomb_event.h"
#include "recomb_prob.h"
#include "MRCA_checks.h"
//

Arg::Arg(int n,double rho,double delta,vector<int> blocks,vector<int> gaps,PopSize * popsize) {
  this->n=n;
  this->rho=rho;
  this->delta=delta;
  this->blocks=blocks;
  this->gaps=gaps;
  this->popsize=popsize;
  changeLT=vector<bool>(blocks.back(),false);
  changeLT[0]=true;
  construct();
}

void Arg::construct() {
  Arg::MRCA M;
  clock_t t1, t2;
  t1=clock();
  s.clear(); //Nodes in the graph
  ages.clear(); //List of ages of all nodes
  clonal.clear(); //Clonal check of all nodes
  int L=blocks.back(); //total number of ancestral sites in the genome
  int b=blocks.size()-1; //Number of blocks

  //Length of genome with gaps inserted between each block
  int G = L;
  for (int i=0;i<b;++i) G += gaps[i];

  const double noStop = 1 - (1/delta); //Geometric rate of elongating recombinant intervals
  const double siteRecomb = rho/(2*G); //Individual site rate of recombinant interval initiation

  int k=n; //Current number of nodes is the number of initial isolates
  vector<int> toCoal;//Contains the list of lines currently in the ARG
  vector<double> recombRates;//Recombination rate of each node
  vector<vector<double> > probStart; //Probability of a recombination interval starting at the start of each interval of ancestral material
  vector<int> blockStarts; //Contains a list of where the blocks start on the genome
  vector<int> blockEnds; //Contains a list of where the blocks end on the genome
  list<list<int> > intervalStarts; //A list of all ancestral interval start sites for each node
  list<list<int> > intervalEnds;//A list of all ancestral interval end sites for each node
  vector<int> totMaterial; //Total amount of ancestral material contained at each node

  //Calculate positions of ancestral blocks with gaps
  blockStarts.push_back(blocks[0]);
  blockEnds.push_back(blocks[1]);
  int gap = 0;
  for (int i=1;i<b;++i){
    gap += gaps[i-1];
    blockStarts.push_back(blocks[i] + gap);
    blockEnds.push_back(blocks[i+1] + gap);
  }
  blockStarts.push_back(G);

  //Set values for the first node to be copied into all other leaves
  toCoal.push_back(0);
  s.push_back(vector<int>(6,-1));
  ages.push_back(0.0);
  clonal.push_back(true);

  //Create ancestral intervals
  intervalStarts.push_back(list<int>());
  intervalEnds.push_back(list<int>());
  for (int i=0;i<b;++i){
    intervalStarts.back().push_back(blockStarts[i]);
    intervalEnds.back().push_back(blockEnds[i]-1);
  }

  //Initialise MRCA struct
  list<int>::iterator itStart = (intervalStarts.back()).begin(), itEnd = (intervalEnds.back()).begin();
  ++itStart;
  //Join together intervals with no gaps
  while (itStart != (intervalStarts.back()).end()){
    if (*itStart <= (*itEnd + 1)){
      itStart = (intervalStarts.back()).erase(itStart);
      itEnd = (intervalEnds.back()).erase(itEnd);
    }else{
      ++itStart;
      ++itEnd;
    }
  }
  itStart = (intervalStarts.back()).begin(), itEnd = (intervalEnds.back()).begin();
  while (itStart != (intervalStarts.back()).end()){
    M.starts.push_back(*itStart);
    M.ends.push_back(*itEnd);
    M.values.push_back(n);//Initialise material for all leaf nodes
    ++itStart;
    ++itEnd;
  }

  //Calculate probability of recombination at each site
  probStart.push_back(vector<double>(b,0.0));
  recombRates.push_back(0.0);
  totMaterial.push_back(0);
  calc_clonalRecomb(G, delta, probStart.back(), recombRates.back(), intervalStarts.back(), intervalEnds.back(), noStop, siteRecomb, totMaterial.back());

  //Set values for initial nodes of the ARG
  for (int i=1;i<n;++i){
    toCoal.push_back(i); //Initial nodes are in the ARG
    s.push_back(vector<int>(6,-1)); //Create new node
    ages.push_back(0.0); //Give all nodes an age of 0
    clonal.push_back(true); //All initial nodes are clonal
    probStart.push_back(vector<double>(b,0.0));
    probStart.back() = probStart[0];
    recombRates.push_back(0.0);
    recombRates.back() = recombRates[0];
    intervalStarts.push_back(list<int>());
    intervalStarts.back() = intervalStarts.front();
    intervalEnds.push_back(list<int>());
    intervalEnds.back() = intervalEnds.front();
    totMaterial.push_back(0);
    totMaterial.back() = totMaterial[0];
  }

  double currentTime=0.0;
  t2=clock();
  cout << "Time spent initialising nodes: " << (double)(t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;
  double split_time=0.0, recomb_probTime=0.0, recomb_intervalTime=0.0, combine_time=0.0, check_time=0.0, remove_time=0.0;

  //Simulate the coalescence-recombination graph
  while (k>1) {
    //Calculate the current rate of recombination
    double currentRecomb = 0.0;
    for (size_t i=0;i<recombRates.size();++i) currentRecomb += recombRates[i];
    //Simulate time to next event exponentially
    currentTime+=gsl_ran_exponential(rng,2.0/(k*(k-1)+(2.0*currentRecomb)));
    //Randomly choose coalescence or recombination

    if (gsl_rng_uniform(rng)<(k*(k-1.0))/(k*(k-1.0)+(2.0*currentRecomb))){
      //Coalescence event
      //Choose two children to coalesce at random
      int i = floor(gsl_rng_uniform(rng)*k);
      int j = i;
      while (j==i) j=floor(gsl_rng_uniform(rng)*k);
      if (i == k-1){
        i=j;
        j=k-1;
      }

      //Create new node
      s.push_back(vector<int>(6,-1));
      //Set two children of new node
      s.back()[0] = toCoal[i];
      s.back()[1] = toCoal[j];
      //Set parent of children
      s[toCoal[i]][2] = s.size()-1;
      s[toCoal[j]][2] = s.size()-1;
      //Add the age of the new node
      ages.push_back(currentTime);
      //Determine whether new node is clonal based on clonal status of either child
      clonal.push_back(clonal[toCoal[i]]||clonal[toCoal[j]]);
      //Add new node to the current ARG
      toCoal[i] = s.size()-1;

      //Test for fully coalesced material
      t1=clock();
      list<list<int> >::iterator itChildStart1 = intervalStarts.begin(), itChildEnd1 = intervalEnds.begin(), itChildStart2 = intervalStarts.begin(), itChildEnd2 = intervalEnds.begin();
      advance(itChildStart1,i);
      advance(itChildEnd1,i);
      advance(itChildStart2,j);
      advance(itChildEnd2,j);
      update_MRCA(M, *itChildStart1, *itChildEnd1, *itChildStart2, *itChildEnd2);
      t2=clock();
      check_time += t2-t1;

      t1=clock();
      combine_ancestries(*itChildStart1, *itChildEnd1, *itChildStart2, *itChildEnd2);
      t2=clock();
      combine_time += t2-t1;

      //Remove the second child from the ARG
      toCoal[j] = toCoal.back();
      toCoal.pop_back();
      recombRates[j] = recombRates.back();
      recombRates.pop_back();
      probStart[j] = probStart.back();
      probStart.pop_back();
      *itChildStart2 = intervalStarts.back();
      intervalStarts.pop_back();
      *itChildEnd2 = intervalEnds.back();
      intervalEnds.pop_back();
      totMaterial[j] = totMaterial.back();
      totMaterial.pop_back();
      --k;

      //Find blocks in MRCA struct that have reached a value of one, fully coalesced
      t1=clock();
      M.itStart = M.starts.begin();
      M.itEnd = M.ends.begin();
      M.itValue = M.values.begin();
      while (M.itValue != (M.values).end()){
        if (*(M.itValue) == 1){
          //MRCA reached for interval, remove this interval from ancestral material
          removeAncMat(*M.itStart, *M.itEnd, *itChildStart1, *itChildEnd1);
          //Remove this interval from the MRCA struct
          M.itStart = (M.starts).erase(M.itStart);
          M.itEnd = (M.ends).erase(M.itEnd);
          M.itValue = (M.values).erase(M.itValue);
          int tempVal = *(M.itValue);
          --M.itValue;
          //Check to see if intervals either side of removed interval can be merged
          if (*(M.itValue) == tempVal){
            M.itStart = (M.starts).erase(M.itStart);
            --M.itEnd;
            M.itEnd = (M.ends).erase(M.itEnd);
            ++M.itEnd;
            M.itValue = (M.values).erase(M.itValue);
            ++M.itValue;
          }else ++M.itValue;
        }else{
          ++M.itStart;
          ++M.itEnd;
          ++M.itValue;
        }
      }
      t2=clock();
      remove_time += (t2-t1);

      t1 = clock();
      //Calculate the recombination rate for the new node
      if (clonal[toCoal[i]] == true){
        calc_clonalRecomb(G, delta, probStart[i], recombRates[i], *itChildStart1, *itChildEnd1, noStop, siteRecomb, totMaterial[i]);
      }else{
        calc_nonClonalRecomb(G, delta, probStart[i], recombRates[i], *itChildStart1, *itChildEnd1, noStop, siteRecomb, totMaterial[i]);
      }
      t2=clock();
      recomb_probTime += (t2-t1);

    }else{
      //Recombination event

      //Choose a child to undergo recombination weighted by its local recombination rate
      double r_1=gsl_rng_uniform(rng);
      int i=0;
      while (r_1>(recombRates[i]/currentRecomb)){
        r_1 -= recombRates[i]/currentRecomb;
        ++i;
      }

      //Choose a start site for recombination based on the probstart vector
      int beg=0, end=0;
      t1=clock();
      list<list<int> >::iterator itChildStart1 = intervalStarts.begin(), itChildEnd1 = intervalEnds.begin();
      advance(itChildStart1,i);
      advance(itChildEnd1,i);
      if (clonal[toCoal[i]] == true){
        choose_clonalRecomb(probStart[i], G, *itChildStart1, *itChildEnd1, beg, end, delta, totMaterial[i], recombRates[i]);
      }else{
        choose_nonClonalRecomb(probStart[i], G, *itChildStart1, *itChildEnd1, beg, end, noStop, totMaterial[i], recombRates[i]);
      }
      t2=clock();
      recomb_intervalTime += (t2-t1);

      //Choose first parent to be clonal if child is clonal
      clonal.push_back(clonal[toCoal[i]]);
      //Other parent is not clonal
      clonal.push_back(false);
      //Check if the local tree changes in this interval
      //Local tree recombination interval relates to the absolute ancestral material without any gaps
      int LTbeg = beg;
      int LTend = end;
      for (int m=0;m<b;++m){
        if ((LTbeg >= blockStarts[m]) && (LTbeg <= blockEnds[m])){
          LTbeg = LTbeg - blockStarts[m] + blocks[m];
          break;
        }
      }
      for (int m=0;m<b;++m){
        if ((LTend >= blockStarts[m]) && (LTend <= blockEnds[m])){
          LTend = LTend - blockStarts[m] + blocks[m];
          break;
        }
      }
      changeLT[LTbeg] = true;
      changeLT[LTend] = true;

      //Add the ages of new recombinant nodes
      ages.push_back(currentTime);
      ages.push_back(currentTime);
      //Add a new node
      s.push_back(vector<int>(6,-1));
      //Add child to new parent
      s.back()[0]=toCoal[i];
      //Add second parent
      s.push_back(vector<int>(6,-1));
      //Add child to new parent
      s.back()[0]=toCoal[i];
      //Add start and end of import
      s.back()[4]=LTbeg;
      s.back()[5]=LTend;
      //Set the parents of the child node
      s[toCoal[i]][2]=s.size()-2;
      s[toCoal[i]][3]=s.size()-1;
      //Put the clonal parent in the ARG
      toCoal[i]=s.size()-2;
      //Add other parent to ARG
      toCoal.push_back(s.size()-1);

      //Set ancestral material of the parents
      t1=clock();
      intervalStarts.push_back(list<int>());
      intervalEnds.push_back(list<int>());
      split_ancestries(*itChildStart1, *itChildEnd1, intervalStarts.back(), intervalEnds.back(), beg, end);
      t2=clock();
      split_time = (t2-t1);

      //Calculate recombination rates and start-point probabilities for the new parents
      t1=clock();
      if (clonal[toCoal[i]] == true){
        calc_clonalRecomb(G, delta, probStart[i], recombRates[i], *itChildStart1, *itChildEnd1, noStop, siteRecomb, totMaterial[i]);
      }else{
        calc_nonClonalRecomb(G, delta, probStart[i], recombRates[i], *itChildStart1, *itChildEnd1, noStop, siteRecomb, totMaterial[i]);
      }
      //And for non-clonal parent
      recombRates.push_back(0.0);
      probStart.push_back(vector<double>());
      totMaterial.push_back(0);
      calc_nonClonalRecomb(G, delta, probStart.back(), recombRates.back(), intervalStarts.back(), intervalEnds.back(), noStop, siteRecomb, totMaterial.back());
      t2=clock();
      recomb_probTime += (t2-t1);
      ++k;
    }
  }
  cout << "Time spent on calculating recombinant interval: " << (double)recomb_intervalTime/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent on splitting ancestries: " << (double)split_time/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent on calculating recomb probabilities: " << (double)recomb_probTime/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent on combining lineages: " << (double)combine_time/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent updating MRCA struct: " << (double)check_time/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent checking MRCA struct: " << (double)remove_time/CLOCKS_PER_SEC << " seconds" << endl;
  if (popsize!=NULL) for (unsigned int i=0;i<ages.size();i++) ages[i]=popsize->convert(ages[i]);
}

Data * Arg::drawData(double theta) {
  string done;
  int L=blocks.back();
  vector<string*> genotypes(s.size(),NULL);
  genotypes[s.size()-1]=new string(L,'N');
  for (int j=0;j<L;++j) genotypes[s.size()-1]->at(j)=floor(gsl_rng_uniform(rng)*4);//Randomly simulate the genome of the MRCA
  for (int i=s.size()-2;i>=0;--i) {
      genotypes[i]=new string(*(genotypes[s[i][2]]));//Copy data from first parent
      //If parent's first child is a leaf, or the genotype of the first parent's first child is NULL AND the same for the first parent's second child:
      if ((s[s[i][2]][0]<0 || genotypes[s[s[i][2]][0]]!=NULL) && (s[s[i][2]][1]<0 || genotypes[s[s[i][2]][1]]!=NULL)) {
          //Delete the genotype of the parent and put the genotype of the parent as done
          delete(genotypes[s[i][2]]);
          genotypes[s[i][2]]=&done;
        }
      if (s[i][3]>0) {//If there is a second parent, copy the imported fragment
          int beg=s[s[i][3]][4];//Start of import
          int  nd=s[s[i][3]][5];//End of import
          if (beg <= nd){
            //Import contained within genome
            for (int j=beg;j<=nd;++j) genotypes[i]->at(j)=genotypes[s[i][3]]->at(j);//Copy the genotypes of the second parent for the imported fragment
          }else{
            //Import wraps around end of genome
            for (int j=0;j<=nd;++j) genotypes[i]->at(j)=genotypes[s[i][3]]->at(j);
            for (int j=beg;j<L;++j) genotypes[i]->at(j)=genotypes[s[i][3]]->at(j);
          }
          if ((s[s[i][3]][0]<0 || genotypes[s[s[i][3]][0]]!=NULL) && (s[s[i][3]][1]<0 || genotypes[s[s[i][3]][1]]!=NULL)) {
            //Second parent's children are either leaves or NULL, therefore finished with parent
              delete(genotypes[s[i][3]]);
              genotypes[s[i][3]]=&done;
            }
        }
      //Add mutations
      int nbmuts=gsl_ran_poisson(rng,theta/2.0*(ages[s[i][2]]-ages[i]));
      for (int m=0;m<nbmuts;m++) {
          int loc=floor(gsl_rng_uniform(rng)*L);
          genotypes[i]->at(loc)=(genotypes[i]->at(loc)+1+(int)floor(gsl_rng_uniform(rng)*3))%4;//Add the mutations randomly based on a poission process
        }
    }
  //Create data object
  Data * data=new Data(n,blocks);
  for (int i=0;i<n;++i) {
      for (int j=0;j<L;++j) data->set_NO_POLY_UPDATE(i,j,genotypes[i]->at(j));
      delete(genotypes[i]);
    }
  return data;
}

string Arg::extractCG() {
  vector<bool> s4(s.size(),false);//Whether to keep a node or not
  vector<vector<int> > s2=s;
  //First add all nodes on clonal branches
  for (unsigned int k=0;k<s4.size();++k) s4[k]=clonal[k];
  //Second remove nodes with a single son, updating the branching matrix accordingly
  for (int i=n;i<(int)s.size();++i) {
      if (s4[i]==false) continue;
      if (s[i][0]<0 || s4[s[i][0]]==false) swap(s[i][0],s[i][1]);//If first child is a leaf or non-clonal swap the nodes
      if (s[i][1]<0 || s4[s[i][1]]==false) {
          s4[i]=false;
          if (s[i][0]>=0) {
              s[s[i][0]][2]=s[i][2];//Set the first parent of the node as the first parent of its first child
              if (s[i][2]>=0) {
                  if (s[s[i][2]][0]==i) s[s[i][2]][0]=s[i][0];//If the node is its own first parent's first child, set the first child of its first parent as its own first child
                  else s[s[i][2]][1]=s[i][0];//Otherwise the second child is its own first child
                }
            }
        }
    }
  //Find new root
  int ii=s.size()-1;while (s4[ii]==false) --ii;//Find the highest clonal node on the ARG
  if (s[ii][0]<0 || s4[s[ii][0]]==false || s[ii][1]<0 || s4[s[ii][1]]==false) s4[ii]=false;//If root's first child is a node OR the first child is non-clonal OR the same for the second child, the root is non-clonal
  while (s4[ii]==false) --ii;//Find new root
  //Construct tree from root
  string str=buildTree(ii).append(";");
  s=s2;
  return str;
}

string Arg::extractLT(int site) {
  vector<bool> s4(s.size(),false);//Whether to keep a node or not
  vector<vector<int> > s2=s;
  //First add all nodes which are ancestral for the given site
  for (int k=0;k<n;++k) s4[k]=true;
  for (unsigned int k=0;k<s4.size()-1;++k) {
      if (s4[k]==false) continue;//Node not included
      if (s[k][3]==-1) s4[s[k][2]]=true;else {//If node has only one child, set the first child as true
        int beg=s[s[k][3]][4];
        int  nd=s[s[k][3]][5];
        //Set equal to true if the site of interest is included in the imported interval
        if (beg <= nd){
          if (site>=beg && site<=nd) s4[s[k][3]]=true;else s4[s[k][2]]=true;
        }else{
          if ((site <= nd) || (site >= beg)) s4[s[k][3]]=true;else s4[s[k][2]]=true;
        }
      }
    }
//Second remove nodes with a single son, updating the branching matrix accordingly
for (int i=n;i<(int)s.size();++i) {
      if (s4[i]==false) continue;//Node not included for site
      if (s[i][2]<0 || s4[s[i][2]]==false) swap(s[i][2],s[i][3]);//If first parent is a leaf or not included, swap parents
      if (s[i][0]<0 || s4[s[i][0]]==false) swap(s[i][0],s[i][1]);//If first child is a leaf or not include, swap children
      if (s[i][1]<0 || s4[s[i][1]]==false) {//If second child is a leaf or not included
          s4[i]=false;//Set current node to false
          if (s[i][0]>=0) {//If there is a first child now
              s[s[i][0]][2]=s[i][2];//Set first parent of child as first parent of current node
              if (s[i][2]>=0) {//If first parent exists
                  if (s[s[i][2]][0]==i) s[s[i][2]][0]=s[i][0];//If the current node is its first parent's first child, set first child as parent's first child
                  else s[s[i][2]][1]=s[i][0];//Otherwise set first child as first parent's second child
                }
            }
        }
    }
  //Find new root
  int ii=s.size()-1;while (s4[ii]==false) --ii;
  if (s[ii][0]<0 || s4[s[ii][0]]==false || s[ii][1]<0 || s4[s[ii][1]]==false) s4[ii]=false;//If children of new root are leaves or not included, remove current root
  while (s4[ii]==false) ii--;//Fidn new root
  //Construct tree from root
  string str=buildTree(ii).append(";");
  s=s2;
  return str;
}

string Arg::buildTree(int r) {
  ostringstream stm;
  if (r<n) {
      stm<<r<<":"<<ages[s[r][2]];
    }else{
      stm<<"("<<buildTree(s[r][0])<<","<<buildTree(s[r][1])<<"):";
      if (s[r][2]<0) stm<<0.0;
      else stm<<ages[s[r][2]]-ages[r];
    }
  return stm.str();
}

void Arg::outputDOT(ostream * out,bool am) {
  *out<<"digraph G {fontsize=5;ranksep=0.02;ratio=fill;size=\"10,10\";"<<endl<<"edge[arrowhead=none];"<<endl;

//Put the sampled individuals on the same rank at the bottom of the graph
  *out<<"{rank=same;";
  for (int i=1;i<=n;i++) *out<<i<<"[shape=point] ";
  *out<<"}"<<endl;
  int lwd=2;//Line width
  vector<vector<bool> > ancmat;
  int L=blocks.back();
  int ma=0;for (unsigned int i=1;i<blocks.size();++i) ma=max(ma,blocks[i]-blocks[i-1]);
  double width=ceil(50.0/ma);//Width of ancestral material segment, either 50 segments per block, or fewer if block length < 50
  int skip=max(1.0,floor(ma/50.0));//Length of each segment in nucleotides
  for (int i=0;i<n;i++) ancmat.push_back(vector<bool>(L,true));
  for (int i=0;i<(int)ages.size();++i) {
      if (!am) {
        //Ancestral material not to be included in the DOT file
          if (clonal[i]) *out<<i+1<<"[shape=point,width=0.00,height=0.00]"<<endl;//Black node if clonal
          else *out<<i+1<<"[shape=point,width=0.00,height=0.00,color=gray]"<<endl;//Gray node if non-clonal
          continue;
        }
      if (i>=n) {
        //Draw the node with the ancestral material included
          ancmat.push_back(ancmat[s[i][0]]);
          if (s[i][1]>=0) {//Father of two sons
              for (int k=0;k<L;++k) ancmat.back()[k]=ancmat[s[i][0]][k]||ancmat[s[i][1]][k];
            } else if (s[s[i][0]][2]==i) //recipient
          for (int k=0;k<L;++k) {
            if (s[s[s[i][0]][3]][4]<=s[s[s[i][0]][3]][5]){
              if (k>s[s[s[i][0]][3]][4]&&k<s[s[s[i][0]][3]][5]) ancmat.back()[k]=false;
            }else{
              if (k<s[s[s[i][0]][3]][5]||k>s[s[s[i][0]][3]][4]) ancmat.back()[k]=false;
            }
          }
          else//donor
          if (s[i][4] <= s[i][5]){
            for (int k=0;k<L;++k) {if (k<=s[i][4]||k>=s[i][5]) ancmat.back()[k]=false;}
          }else{
            for (int k=0;k<L;++k) {if (k>=s[i][5]&&k<=s[i][4]) ancmat.back()[k]=false;}
          }
        }
        //Draw each node as bars with ancestral material colourd gray on white background
      *out<<i+1<<"[shape=plaintext,label=<<table CELLBORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"0\" BORDER=\"0\">";
      //Draw white border at top of box
      *out<<"<tr><td HEIGHT=\""<<lwd<<"\" COLSPAN=\""<<floor(ma/skip)+4<<"\" bgcolor=\"white\"></td></tr><tr>";
      for (unsigned int j=1;j<blocks.size();++j) {
        //Draw black line at top of box
          *out<<"<td HEIGHT=\""<<lwd<<"\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td><td bgcolor=\"black\" HEIGHT=\""<<lwd<<"\" COLSPAN=\""<<floor((blocks[j]-blocks[j-1])/skip)+2<<"\"></td>";
      }
      *out <<"<td HEIGHT=\""<<lwd<<"\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td></tr><tr>";
      for (unsigned int j=1;j<blocks.size();++j){
        //Draw vertical line at left of box
          *out<<"<td HEIGHT=\"10\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td><td bgcolor=\"black\" WIDTH=\""<<lwd<<"\" HEIGHT=\"10\"></td>";
        //Draw vertical grey or white lines to represent ancestral material or lack of
          for (int k=blocks[j-1];k<blocks[j];k+=skip)
            if (ancmat[i][k])
              *out<<"<td bgcolor=\"grey\" HEIGHT=\"10\" WIDTH=\""<<width<<"\"></td>";
            else
              *out<<"<td HEIGHT=\"10\" WIDTH=\""<<width<<"\" bgcolor=\"white\"></td>";
        //Draw vertical black line at left hand end of box
          *out<<"<td bgcolor=\"black\" WIDTH=\""<<lwd<<"\" HEIGHT=\"10\"></td>";
      }
        //Draw white border at end
          *out<<"<td bgcolor=\"white\" WIDTH=\""<<lwd<<"\" HEIGHT=\"10\"></td></tr><tr>";
      for (unsigned int j=1;j<blocks.size();++j){
        //Draw black line along bottom of box
          *out<<"<td HEIGHT=\""<<lwd<<"\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td><td bgcolor=\"black\" HEIGHT=\""<<lwd<<"\" COLSPAN=\""<<floor((blocks[j]-blocks[j-1])/skip)+2<<"\"></td>";
      }
      *out<<"<td HEIGHT=\""<<lwd<<"\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td></tr>";
        //Draw white border at bottom of box
      *out<<"<tr><td HEIGHT=\"5\" COLSPAN=\""<<floor(ma/skip)+4<<"\" bgcolor=\"white\"></td></tr>";
      *out<<"</table>>]"<<endl;
    }

//For standard colouring of the edges
  for (unsigned int a=0;a<s.size();++a) {
      //Clonal frame edges
      if (s[a][2]>=0) {
          int j=s[a][2];
          //if ages(a)==-1;continue;end
          if (clonal[j] && clonal[a])
            *out<<j+1<<" -> "<<a+1<<"[style=bold]"<<endl;
          else
            *out<<j+1<<" -> "<<a+1<<"[color=gray]"<<endl;
        }

      //Recombination edges
      if (s[a][3]>=0) {
          int j=s[a][3];
          //int o=s[a][2];
          //if ages(o)==-1;continue;end
          if (clonal[j] && clonal[a])
            *out<<j+1<<" -> "<<a+1<<"[style=bold]"<<endl;
          else
            *out<<j+1<<" -> "<<a+1<<"[color=gray]"<<endl;
        }
    }
  *out<<"}"<<endl;
}

void Arg::outputLOCAL(ostream * out) {
  unsigned int i=0;
  while (1) {
      string tree=extractLT(i);
      int n=0;
      while (i+1<changeLT.size() && changeLT[i+1]==false) {++i;++n;}
      *out<<"["<<n+1<<"]"<<tree<<endl;
      ++i;
      if (i==changeLT.size()) break;
    }
}
