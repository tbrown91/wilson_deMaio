#include "arg.h"
#include <math.h>
#include <ctime>
#include "coal_event.h"
#include "recomb_event.h"
#include "recomb_prob.h"
#include "calc_intervals.h"
//#include "updateMRCA.h"
//
struct MRCA {
    list<int> starts;
    list<int> ends;
    list<int> values;
    list<int>::iterator itStart;
    list<int>::iterator itEnd;
    list<int>::iterator itValue;
};

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

void removeAncMat(const int start, const int end, vector<int> &starts, vector<int> &ends){
  //Reove MRCA material from chosen nodes
  //cout << "Interval to remove: " << start << " " << end << endl;
  vector<int> tempStarts;
  vector<int> tempEnds;
  for (int a=0;a<int(starts.size());++a){
    if ((end < starts[a]) || (start > ends[a])){
      tempStarts.push_back(starts[a]);
      tempEnds.push_back(ends[a]);
    }else if ((start <= starts[a]) && (end >= starts[a])){
      if (end < ends[a]){
        tempStarts.push_back(end+1);
        tempEnds.push_back(ends[a]);
      }
      //Otherwise remove interval from ancestral material
    }else if ((end >= ends[a]) && (start <= ends[a])){
      tempStarts.push_back(starts[a]);
      tempEnds.push_back(start-1);
    }else if ((start > starts[a]) && (end < ends[a])){
      tempStarts.push_back(starts[a]);
      tempEnds.push_back(start-1);
      tempStarts.push_back(end+1);
      tempEnds.push_back(ends[a]);
    }
  }
  starts = tempStarts;
  ends = tempEnds;
}

void modifyMRCA(MRCA& M, const int start, const int end){
//update the list of MRCA intervals by subtracting an interval overlapping in the two coalesceing lineages
  //cout << "Overlap interval:" << endl;
  //cout << start << " " << end << endl;
  M.itStart=M.starts.begin();
  M.itEnd=M.ends.begin();
  M.itValue=M.values.begin();
	if ((M).itStart==(M).starts.end()){
		cout << "Error : this part should be in the MRCA vector, but it isn't (1) " << *M.itStart << " " << *M.starts.end() << endl;
	 	return;
   //exit();
	}
	while (*((M).itEnd) < start) {
    //cout << "Start: " << *M.itStart << " End: " << *M.itEnd << " Value: " << *M.itValue << endl;
		++((M).itStart);
		++((M).itEnd);
		++((M).itValue);
		if ((M).itStart==((M).starts).end()){
			cout << "Error : this part should be in the MRCA vector, but it isn't (2)" << *M.itStart << " " << *M.starts.end() << endl;
    	//exit();
      break;
		}
	}
  //cout << "Start: " << *M.itStart << " End: " << *M.itEnd << " Value: " << *M.itValue << endl;

	if (start<*((M).itStart)){
		cout << "Error : this part should be in the MRCA vector, but it isn't (3)" << start << " " << *M.itStart << endl;
  	//exit();
    return;
	}

	//first element in MRCA list that overlaps, but starting values do not coincide
	if (((M).itStart!=(M).starts.end()) && *((M).itStart)<start){
    ++((M).itStart);
    (M).starts.insert((M).itStart,start);
    --((M).itStart);
		(M).ends.insert((M).itEnd,start-1);
		(M).values.insert((M).itValue,*(M).itValue);
	}

  	//iteratively look at all intervals in the MRCA structure overlapping completely the given interval
	while ((((M).itStart)!=(M).starts.end()) && (*((M).itEnd)<=end)) {
    --*((M).itValue);
		++((M).itStart);
		++((M).itEnd);
		++((M).itValue);
	}

  // look at the final interval in the MRCA, if present, that overlaps the given interval but has end after the end of the given interval
  if (((M).itStart!=(M).starts.end()) && (*((M).itStart)<=end)){
    ++((M).itStart);
    (M).starts.insert((M).itStart,end+1);
    (M).ends.insert((M).itEnd,end);
    (M).values.insert((M).itValue,*((M).itValue));
    --(M).itValue;
    --*((M).itValue);
  }
}

void Arg::construct() {
  MRCA M;
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

  const double noStop = 1 - (1/delta);
  const double siteRecomb = rho/(2*G);

  int k=n; //Current number of nodes is the number of initial isolates
  vector<int> toCoal;//Contains the list of lines currently in the ARG
  vector<vector<bool> > toCoalAncMat;//Describes the ancestral material of these lines
  vector<double> recombRates;//<Recombination rate of each node
  vector<vector<double> > probStart; //Probability of a recombination interval starting at each site of the genome
  vector<int> blockStarts; //Contains a list of where the blocks start on the genome
  vector<int> blockEnds; //Contains a list of where the blocks end on the genome
  vector<vector<int> > intervalStarts; //A list of all ancestral interval start sites for each node
  vector<vector<int> > intervalEnds;//A list of all ancestral interval end sites for each node

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
  s.push_back(vector<int>(6.-1));
  ages.push_back(0.0);
  clonal.push_back(true);
  //Create ancestral intervals
  intervalStarts.push_back(vector<int>(NULL));
  intervalEnds.push_back(vector<int>(NULL));
  for (int i=0;i<b;++i){
    intervalStarts.back().push_back(blockStarts[i]);
    intervalEnds.back().push_back(blockEnds[i]);
  }
  //Initialise MRCA struct
  for (int i=0;i<b;++i){
    M.starts.push_back(intervalStarts.back()[i]);
    M.ends.push_back(intervalEnds.back()[i]);
    M.values.push_back(n);//Initialise material for all leaf nodes
  }
  //Calculate probability of recombination at each site
  probStart.push_back(vector<double>(G,0.0));
  recombRates.push_back(0.0);
  calc_clonalRecomb(G, delta, probStart.back(), recombRates.back(), intervalStarts.back(), intervalEnds.back(), noStop, siteRecomb);

  //Set values for initial nodes of the ARG
  for (int i=1;i<n;i++) {
      toCoal.push_back(i); //Initial nodes are in the ARG
      s.push_back(vector<int>(6,-1)); //Create new node
      ages.push_back(0.0); //Give all nodes an age of 0
      clonal.push_back(true); //All initial nodes are clonal
      //Set the probability of starting recombination at each site of the genome
      probStart.push_back(vector<double>(G,0.0));
      probStart.back() = probStart[0];
      recombRates.push_back(0.0);
      recombRates.back() = recombRates[0];
      intervalStarts.push_back(vector<int>(NULL));
      intervalStarts.back() = intervalStarts[0];
      intervalEnds.push_back(vector<int>(NULL));
      intervalEnds.back() = intervalEnds[0];
  }

  double currentTime=0.0;
  //Simulate the coalescence-recombination graph
  t2=clock();
  cout << "Time spent initialising nodes: " << (double)(t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;
  double split_time=0.0, recomb_probTime=0.0, recomb_intervalTime=0.0, combine_time=0.0, remove_time=0.0;
  while (k>1) {
    //cout << "Nodes remaining: " << k << endl;
    //Calculate the current rate of recombination
    double currentRecomb = 0.0;
    for (int i=0;i<int(recombRates.size());++i) currentRecomb += recombRates[i];
    //Simulate time to next event exponentially
    currentTime+=gsl_ran_exponential(rng,2.0/(k*(k-1)+(2.0*currentRecomb)));
    //Randomly choose coalescence or recombination
    if (gsl_rng_uniform(rng)<1.0*(k*(k-1))/(k*(k-1)+(2.0*currentRecomb))){
      //cout << "Coal" << endl;
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
      {
        int index1=0, index2=0;
        int b1 = intervalStarts[i].size(), b2 = intervalStarts[j].size();
        // cout << "Children:" << endl;
        // for (int a=0;a<b1;++a) cout << intervalStarts[i][a] << " " << intervalEnds[i][a] << " ";
        // cout << endl;
        // for (int a=0;a<b2;++a) cout << intervalStarts[j][a] << " " << intervalEnds[j][a] << " ";
        // cout << endl;
        int currentStart1=0, currentStart2=0;
        if ((index1 != b1) && (index2 != b2)) currentStart1=intervalStarts[i][index1], currentStart2=intervalStarts[j][index2];
        while ((index1 != b1) && (index2 != b2)){
          if ((currentStart1 < currentStart2) && (intervalEnds[i][index1]>=currentStart2)){ //overlap
            currentStart1=currentStart2;
          }else if ((currentStart2 < currentStart1) && (intervalEnds[j][index2]>=currentStart1)){ //overlap
            currentStart2=currentStart1;
          }else if ((currentStart1 < currentStart2) && (intervalEnds[i][index1]<currentStart2)){ //no overlap
            ++index1;
            if (index1<b1) currentStart1=intervalStarts[i][index1];
          }else if ((currentStart2 < currentStart1) && (intervalEnds[j][index2]<currentStart1)){ //no overlap
            ++index2;
            if (index2<b2) currentStart2=intervalStarts[j][index2];
          }else if ((currentStart1 == currentStart2) && (intervalEnds[i][index1]<intervalEnds[j][index2])){ //overlap, same start
            modifyMRCA(M, currentStart1, intervalEnds[i][index1]);
            currentStart2=intervalEnds[i][index1]+1;
            ++index1;
            if (index1<b1) currentStart1=intervalStarts[i][index1];
          }else if ((currentStart1 == currentStart2) && (intervalEnds[j][index2]<intervalEnds[i][index1])){ //overlap, same start
            modifyMRCA(M, currentStart1, intervalEnds[j][index2]);
            currentStart1=intervalEnds[j][index2]+1;
            ++index2;
            if (index2<b2) currentStart2=intervalStarts[j][index2];
          }else if ((currentStart1 == currentStart2) && (intervalEnds[j][index2]==intervalEnds[i][index1])){ //overlap, same start, same end
            modifyMRCA(M, currentStart1, intervalEnds[i][index1]);
            ++index2;
            if (index2<b2) currentStart2=intervalStarts[j][index2];
            ++index1;
            if (index1<b1) currentStart1=intervalStarts[i][index1];
          }else{
            cout << "Error : this case should not happen in updating MRCA vector" << endl;
            //exit();
            break;
          }
        }
      }
      //Find blocks in MRCA struct that have reached a value of zero
      M.itStart = M.starts.begin();
      M.itEnd = M.ends.begin();
      M.itValue = M.values.begin();
      while (M.itStart != M.starts.end()){
        if (*M.itValue == 1){
          removeAncMat(*M.itStart, *M.itEnd, intervalStarts[i], intervalEnds[i]);
          removeAncMat(*M.itStart, *M.itEnd, intervalStarts[j], intervalEnds[j]);
          *M.itValue = -1;
        }
        ++M.itStart;
        ++M.itEnd;
        ++M.itValue;
      }
      t2=clock();
      remove_time += t2-t1;
      t1=clock();
      combine_ancestries(intervalStarts[i], intervalEnds[i], intervalStarts[j], intervalEnds[j]);
      t2=clock();
      combine_time += t2-t1;
      //If material has been removed from the lineage, recalculate the intervals
      t1 = clock();
      //Calculate the recombination rate for the new node
      if (clonal[toCoal[i]] == true){
        calc_clonalRecomb(G, delta, probStart[i], recombRates[i], intervalStarts[i], intervalEnds[i], noStop, siteRecomb);
      }else{
        calc_nonClonalRecomb(G, delta, probStart[i], recombRates[i], intervalStarts[i], intervalEnds[i], noStop, siteRecomb);
      }
      t2=clock();
      recomb_probTime += (t2-t1);
      //Remove the second child from the ARG
      toCoal[j] = toCoal.back();
      toCoal.pop_back();
      recombRates[j] = recombRates.back();
      recombRates.pop_back();
      probStart[j] = probStart.back();
      probStart.pop_back();
      intervalStarts[j] = intervalStarts.back();
      intervalStarts.pop_back();
      intervalEnds[j] = intervalEnds.back();
      intervalEnds.pop_back();
      --k;
    }else{
      //cout << "Recomb" << endl;
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
      if (clonal[toCoal[i]] == true){
        choose_clonalRecomb(probStart[i], G, intervalStarts[i], intervalEnds[i], beg, end, delta);
      }else{
        choose_nonClonalRecomb(probStart[i], G, intervalStarts[i], intervalEnds[i], beg, end, noStop);
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
        if ((LTbeg >= blockStarts[m]) && (LTbeg < blockEnds[m])){
          LTbeg = LTbeg - blockStarts[m] + blocks[m];
          break;
        }
      }
      for (int m=0;m<b;++m){
        if ((LTend >= blockStarts[m]) && (LTend < blockStarts[m+1])){
          if (LTend >= blockEnds[m]) LTend = blocks[m+1]-1;
          else LTend = LTend - blockStarts[m] + blocks[m];
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
      intervalStarts.push_back(vector<int>(NULL));
      intervalEnds.push_back(vector<int>(NULL));
      split_ancestries(intervalStarts[i],intervalEnds[i],intervalStarts.back(), intervalEnds.back(), beg, end);
      t2=clock();
      split_time = (t2-t1);
      //Calculate recombination rates and start-point probabilities for the new parents
      t1=clock();
      if (clonal[toCoal[i]] == true){
        calc_clonalRecomb(G, delta, probStart[i], recombRates[i], intervalStarts[i], intervalEnds[i], noStop, siteRecomb);
      }else{
        calc_nonClonalRecomb(G, delta, probStart[i], recombRates[i], intervalStarts[i], intervalEnds[i], noStop, siteRecomb);
      }
      //And for non-clonal parent
      recombRates.push_back(0.0);
      probStart.push_back(vector<double>(G,0.0));
      calc_nonClonalRecomb(G, delta, probStart.back(), recombRates.back(), intervalStarts.back(), intervalEnds.back(), noStop, siteRecomb);
      t2=clock();
      recomb_probTime += (t2-t1);
      ++k;
    }
  }
  cout << "Time spent on calculating recombinant interval: " << (double)recomb_intervalTime/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent on splitting ancestries: " << (double)split_time/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent on calculating recomb probabilities: " << (double)recomb_probTime/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent on combining lineages: " << (double)combine_time/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent checking for fully coalesced material: " << (double)remove_time/CLOCKS_PER_SEC << " seconds" << endl;
  if (popsize!=NULL) for (unsigned int i=0;i<ages.size();i++) ages[i]=popsize->convert(ages[i]);
}

Data * Arg::drawData(double theta) {
  string done;
  int L=blocks.back();
  vector<string*> genotypes(s.size(),NULL);
  genotypes[s.size()-1]=new string(L,'N');
  for (int j=0;j<L;++j) genotypes[s.size()-1]->at(j)=floor(gsl_rng_uniform(rng)*4);
  for (int i=s.size()-2;i>=0;--i) {
      genotypes[i]=new string(*(genotypes[s[i][2]]));//Copy data from first parent
      if ((s[s[i][2]][0]<0 || genotypes[s[s[i][2]][0]]!=NULL) && (s[s[i][2]][1]<0 || genotypes[s[s[i][2]][1]]!=NULL)) {
          delete(genotypes[s[i][2]]);
          genotypes[s[i][2]]=&done;
        }
      if (s[i][3]>0) {//If there is a second parent, copy the imported fragment
          int beg=s[s[i][3]][4];
          int  nd=s[s[i][3]][5];
          for (int j=beg;j<nd;j++) genotypes[i]->at(j)=genotypes[s[i][3]]->at(j);
          if ((s[s[i][3]][0]<0 || genotypes[s[s[i][3]][0]]!=NULL) && (s[s[i][3]][1]<0 || genotypes[s[s[i][3]][1]]!=NULL)) {
              delete(genotypes[s[i][3]]);
              genotypes[s[i][3]]=&done;
            }
        }
      //Add mutations
      int nbmuts=gsl_ran_poisson(rng,theta/2.0*(ages[s[i][2]]-ages[i]));
      for (int m=0;m<nbmuts;m++) {
          int loc=floor(gsl_rng_uniform(rng)*L);
          genotypes[i]->at(loc)=(genotypes[i]->at(loc)+1+(int)floor(gsl_rng_uniform(rng)*3))%4;
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
      if (s[i][0]<0 || s4[s[i][0]]==false) swap(s[i][0],s[i][1]);
      if (s[i][1]<0 || s4[s[i][1]]==false) {
          s4[i]=false;
          if (s[i][0]>=0) {
              s[s[i][0]][2]=s[i][2];
              if (s[i][2]>=0) {
                  if (s[s[i][2]][0]==i) s[s[i][2]][0]=s[i][0];
                  else s[s[i][2]][1]=s[i][0];
                }
            }
        }
    }
  //Find new root
  int ii=s.size()-1;while (s4[ii]==false) --ii;
  if (s[ii][0]<0 || s4[s[ii][0]]==false || s[ii][1]<0 || s4[s[ii][1]]==false) s4[ii]=false;
  while (s4[ii]==false) --ii;
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
      if (s4[k]==false) continue;
      if (s[k][3]==-1) s4[s[k][2]]=true;else {
        int beg=s[s[k][3]][4];
        int  nd=s[s[k][3]][5];
        if (site>=beg && site<nd) s4[s[k][3]]=true;else s4[s[k][2]]=true;
      }
    }
//Second remove nodes with a single son, updating the branching matrix accordingly
for (int i=n;i<(int)s.size();++i) {
      if (s4[i]==false) continue;
      if (s[i][2]<0 || s4[s[i][2]]==false) swap(s[i][2],s[i][3]);
      if (s[i][0]<0 || s4[s[i][0]]==false) swap(s[i][0],s[i][1]);
      if (s[i][1]<0 || s4[s[i][1]]==false) {
          s4[i]=false;
          if (s[i][0]>=0) {
              s[s[i][0]][2]=s[i][2];
              if (s[i][2]>=0) {
                  if (s[s[i][2]][0]==i) s[s[i][2]][0]=s[i][0];
                  else s[s[i][2]][1]=s[i][0];
                }
            }
        }
    }
  //Find new root
  int ii=s.size()-1;while (s4[ii]==false) --ii;
  if (s[ii][0]<0 || s4[s[ii][0]]==false || s[ii][1]<0 || s4[s[ii][1]]==false) s4[ii]=false;
  while (s4[ii]==false) ii--;
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
  int lwd=2;
  vector<vector<bool> > ancmat;
  int L=blocks.back();
  int ma=0;for (unsigned int i=1;i<blocks.size();++i) ma=max(ma,blocks[i]-blocks[i-1]);
  double width=ceil(50.0/ma);
  int skip=max(1.0,floor(ma/50.0));
  for (int i=0;i<n;i++) ancmat.push_back(vector<bool>(L,true));
  for (int i=0;i<(int)ages.size();++i) {
      if (!am) {
          if (clonal[i]) *out<<i+1<<"[shape=point,width=0.00,height=0.00]"<<endl;
          else *out<<i+1<<"[shape=point,width=0.00,height=0.00,color=gray]"<<endl;
          continue;
        }
      if (i>=n) {
          ancmat.push_back(ancmat[s[i][0]]);
          if (s[i][1]>=0) {//Father of two sons
              for (int k=0;k<L;++k) ancmat.back()[k]=ancmat[s[i][0]][k]||ancmat[s[i][1]][k];
            } else if (s[s[i][0]][2]==i) //recipient
          for (int k=0;k<L;++k) {if (k>=s[s[s[i][0]][3]][4]&&k<s[s[s[i][0]][3]][5]) ancmat.back()[k]=false;}
          else//donor
          for (int k=0;k<L;++k) {if (k<s[i][4]||k>=s[i][5]) ancmat.back()[k]=false;}
        }
      *out<<i+1<<"[shape=plaintext,label=<<table CELLBORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"0\" BORDER=\"0\">";
      *out<<"<tr><td HEIGHT=\""<<lwd<<"\" COLSPAN=\""<<floor(ma/skip)+4<<"\" bgcolor=\"white\"></td></tr>";
      for (unsigned int j=1;j<blocks.size();++j) {
          *out<<"<tr><td HEIGHT=\""<<lwd<<"\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td><td bgcolor=\"black\" HEIGHT=\""<<lwd<<"\" COLSPAN=\""<<floor((blocks[j]-blocks[j-1])/skip)+2<<"\"></td><td HEIGHT=\""<<lwd<<"\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td></tr>";
          *out<<"<tr><td HEIGHT=\"10\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td><td bgcolor=\"black\" WIDTH=\""<<lwd<<"\" HEIGHT=\"10\"></td>";
          for (int k=blocks[j-1];k<blocks[j];k+=skip)
            if (ancmat[i][k])
              *out<<"<td bgcolor=\"grey\" HEIGHT=\"10\" WIDTH=\""<<width<<"\"></td>";
            else
              *out<<"<td HEIGHT=\"10\" WIDTH=\""<<width<<"\" bgcolor=\"white\"></td>";
          *out<<"<td bgcolor=\"black\" WIDTH=\""<<lwd<<"\" HEIGHT=\"10\"></td>";
          *out<<"<td bgcolor=\"white\" WIDTH=\""<<lwd<<"\" HEIGHT=\"10\"></td></tr>";
          *out<<"<tr><td HEIGHT=\""<<lwd<<"\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td><td bgcolor=\"black\" HEIGHT=\""<<lwd<<"\" COLSPAN=\""<<floor((blocks[j]-blocks[j-1])/skip)+2<<"\"></td><td HEIGHT=\""<<lwd<<"\" WIDTH=\""<<lwd<<"\" bgcolor=\"white\"></td></tr>";
          *out<<"<tr><td HEIGHT=\"5\" COLSPAN=\""<<floor(ma/skip)+4<<"\" bgcolor=\"white\"></td></tr>";
        }
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
