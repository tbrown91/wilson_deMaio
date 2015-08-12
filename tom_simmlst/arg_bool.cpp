#include "arg.h"
#include <math.h>
#include <ctime>
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

void calc_intervals(const vector<bool> &AncMat, vector<int> &intervalStarts, vector<int> &intervalEnds, const vector<int> &blockStarts, const vector<int> &blockEnds, const int b){
  //Calculate the intervals of ancestral material for the given node
  if (b == 1){
    //Only one block of ancestral material
    if (AncMat[blockStarts[0]] == true){
      intervalStarts.push_back(blockStarts[0]);
      //Check for interval of length 1
      if (AncMat[blockStarts[0]+1] == false) intervalEnds.push_back(blockStarts[0]);
    }
    for (int a=blockStarts[0]+1;a<blockEnds[0]-1;++a){
      if ((AncMat[a] == true) && (AncMat[a-1] == false)){
        intervalStarts.push_back(a);
        //Check for interval of length 1
        if (AncMat[a+1] == false) intervalEnds.push_back(a);
      }else if ((AncMat[a] == true) && (AncMat[a+1] == false)) intervalEnds.push_back(a);
    }
    //Check last element
    if (AncMat[blockEnds[0]-1] == true){
      intervalEnds.push_back(blockEnds[0]-1);
      //Check for interval of length 1
      if (AncMat[blockEnds[0]-2] == false) intervalStarts.push_back(blockEnds[0]-1);
    }
  }else{
    //More than one ancestral block
    //First ancestral block
    if (AncMat[blockStarts[0]] == true){
      intervalStarts.push_back(blockStarts[blockStarts[0]]);
      //Check for interval of length 1
      if (AncMat[blockStarts[0]+1] == false) intervalEnds.push_back(blockStarts[0]);
    }
    for (int a=blockStarts[0]+1;a<blockEnds[0];++a){
      if ((AncMat[a] == true) && (AncMat[a-1] == false)){
        intervalStarts.push_back(a);
        //Check for interval of length 1
        if (AncMat[a+1] == false) intervalEnds.push_back(a);
      }else if ((AncMat[a] == true) && (AncMat[a+1] == false)) intervalEnds.push_back(a);
    }

    //Middle ancestral blocks
    for (int m=1;m<b-1;++m){
      //Find start and end points
      for (int a=blockStarts[m];a<blockEnds[m];++a){
        if ((AncMat[a] == true) && (AncMat[a-1] == false)){
          intervalStarts.push_back(a);
          //Check for interval of length 1
          if (AncMat[a+1] == false) intervalEnds.push_back(a);
        }else if ((AncMat[a] == true) && (AncMat[a+1] == false)) intervalEnds.push_back(a);
      }
    }
    //Last ancestral block
    for (int a=blockStarts[b-1];a<blockEnds[b-1]-1;++a){
      if ((AncMat[a] == true) && (AncMat[a-1] == false)){
        intervalStarts.push_back(a);
        //Check for interval of length 1
        if (AncMat[a+1] == false) intervalEnds.push_back(a);
      }else if ((AncMat[a] == true) && (AncMat[a+1] == false)) intervalEnds.push_back(a);
    }
    //Check last element of block
    if (AncMat[blockEnds[b-1]-1] == true){
      intervalEnds.push_back(blockEnds[b-1]-1);
      //Check for interval of length 1
      if (AncMat[blockEnds[b-1]-2] == false) intervalStarts.push_back(blockEnds[b-1]-1);
    }
  }
}

void calc_nonClonalRecomb(const int G, const double delta, vector<double> &prob, double &recombRate, const vector<int> &starts, const vector<int> &ends, const double noStop, const double siteRecomb){
  //Calculate the recombination rate and probability of recombination beginning at each site for a non-clonal frame lineage
  //Total amount of ancestral material in the node
  int totMaterial = 0;
  int b = starts.size();
  //Reset recombination rate and probability of recombination along the genome
  prob = vector<double>(G,0.0);
  recombRate = 0;
  for (int m=0;m<b;++m){
    totMaterial += ends[m] - starts[m] + 1;
  }
  if (totMaterial > 1){
    if ((starts[0] == 0) && (ends.back() == G-1)){
      //Last interval wraps around end of genome
      //Ignore first interval start site
      for (int ii=1;ii<b;++ii){
        //Calculate recombination rate at every other start site
        prob[starts[ii]] = siteRecomb*delta*((1-pow(noStop,starts[ii]-ends[ii-1]))*(1-pow(noStop,G+ends[ii-1]-starts[ii])));
        recombRate += prob[starts[ii]];
      }
      //Calculate recombination rate of all sites in intervals past the first site
      recombRate += siteRecomb*(totMaterial-(b-1))*(1-pow(noStop,G-1));
      double recombProb = siteRecomb*(1-pow(noStop,G-1))/recombRate;
      //Calculate probability for first interval
      for (int a=starts[0];a<=ends[0];++a) prob[a] = recombProb;
      for (int ii=1;ii<b;++ii){
        prob[starts[ii]] /= recombRate;
        for (int a=starts[ii]+1;a<=ends[ii];++a) prob[a] = recombProb;
      }
    }else{
      //All intervals are contained within the indices 0 and G-1
      //Calculate recombination rate at first intervals start site
      prob[starts[0]] = siteRecomb*delta*(1-pow(noStop,G+starts[0]-ends.back()))*(1-pow(noStop,ends.back()-starts[0]));
      recombRate += prob[starts[0]];
      for (int ii=1;ii<b;++ii){
        //Calculate recombination rate at every other start site
        prob[starts[ii]] = siteRecomb*delta*((1-pow(noStop,starts[ii]-ends[ii-1]))*(1-pow(noStop,G+ends[ii-1]-starts[ii])));
        recombRate += prob[starts[ii]];
      }
      //Calculate recombination rate of all sites in intervals past the first site
      recombRate += siteRecomb*(totMaterial-b)*(1-pow(noStop,G-1));
      //Normalise the probabilities to the total recombination rate
      double recombProb = siteRecomb*(1-pow(noStop,G-1))/recombRate;
      for (int ii=0;ii<b;++ii){
        prob[starts[ii]] /= recombRate;
        for (int a=starts[ii]+1;a<=ends[ii];++a) prob[a] = recombProb;
      }
    }
  }
}

void calc_clonalRecomb(const int G, const double delta, vector<double> &prob, double &recombRate, const vector<int> &starts, const vector<int> &ends, const double noStop, const double siteRecomb){
  //Total amount of ancestral material in the node
  int totMaterial = 0;
  int b = starts.size();
  //Reset recombination rate and probability of recombination along the genome
  prob = vector<double>(G,0.0);
  recombRate = 0;
  for (int m=0;m<b;++m){
    totMaterial += ends[m] - starts[m] + 1;
  }
  if (totMaterial > 1){
    if ((starts[0] == 0) && (ends.back() == G-1)){
      //Last interval wraps around end of genome
      //Ignore first interval start site
      for (int ii=1;ii<b;++ii){
        //Calculate recombination rate at every other start site
        prob[starts[ii]] = siteRecomb*delta*(1-pow(noStop,starts[ii]-ends[ii-1]));
        recombRate += prob[starts[ii]];
      }
      //Calculate recombination rate of all sites in intervals past the first site
      recombRate += siteRecomb*(totMaterial-(b-1));
      double recombProb = siteRecomb/recombRate;
      //Calculate probability for first interval
      for (int a=starts[0];a<=ends[0];++a) prob[a] = recombProb;
      for (int ii=1;ii<b;++ii){
        prob[starts[ii]] /= recombRate;
        for (int a=starts[ii]+1;a<=ends[ii];++a) prob[a] = recombProb;
      }
    }else{
      //All intervals are contained within the indices 0 and G-1
      //Calculate recombination rate at first intervals start site
      prob[starts[0]] = siteRecomb*delta*(1-pow(noStop,G+starts[0]-ends.back()));
      recombRate += prob[starts[0]];
      for (int ii=1;ii<b;++ii){
        //Calculate recombination rate at every other start site
        prob[starts[ii]] = siteRecomb*delta*(1-pow(noStop,starts[ii]-ends[ii-1]));
        recombRate += prob[starts[ii]];
      }
      //Calculate recombination rate of all sites in intervals past the first site
      recombRate += siteRecomb*(totMaterial-b);
      //Normalise the probabilities to the total recombination rate
      double recombProb = siteRecomb/recombRate;
      for (int ii=0;ii<b;++ii){
        prob[starts[ii]] /= recombRate;
        for (int a=starts[ii]+1;a<=ends[ii];++a) prob[a] = recombProb;
      }
    }
  }
}

void choose_nonClonalRecomb(const vector<double> &prob, const int G, const vector<int> &starts, const vector<int> &ends, int &beg, int &end, const double noStop){
  //Choose a recombination interval for the chosen lineage which is non-clonal
  int b=starts.size();
  //Choose a start site of recombination at random
  double r_2=gsl_rng_uniform(rng);
  beg = 0;
  while (r_2 > prob[beg]){
    r_2 -= prob[beg];
    ++beg;
  }
  if ((starts[0] == 0) && (ends.back() == G-1)){
    //Last interval wraps around end of genome
    //Check if beginning of recombination is at an ancestral interval start site
    int index = -1;
    for (int m=1;m<b;++m){
      if (beg == starts[m]) index = m;
      break;
    }
    if (index > 0){
      //Recombination interval begins at the start of other interval
      //Simulate the end-point via a truncated geometric distribution
      double r_3 = gsl_rng_uniform(rng);
      int len = floor(log(1-r_3*(1-pow(noStop,G+ends[index-1]-starts[index])))/log(noStop));
      end = (beg + len) % G;
    }else{
      //Recombination interval begins within an ancestral interval
      //Simulate the end-point via a truncated geometric distribution
      double r_3 = gsl_rng_uniform(rng);
      int len = floor(log(1-r_3*(1-pow(noStop,G-1)))/log(noStop));
      end = (beg + len) % G;
    }
    //Check if end of recombinant interval falls between ancestral intervals
    for (int a=1;a<b;++a){
      if ((end >= ends[a-1]) && (end < starts[a])){
        end == ends[a-1];
        break;
      }
    }
  }else{
    //All intervals are contained within the "linear" genome
    //Check if beginning of recombination is at an ancestral interval start site
    int index = -1;
    for (int m=0;m<b;++m){
      if (beg == starts[m]) index = m;
      break;
    }
    if (index == 0){
      //Recombination interval begins at the start of first ancestral interval
      //Simulate the end-point via a truncated geometric distribution
      double r_3 = gsl_rng_uniform(rng);
      int len = floor(log(1-r_3*(1-pow(noStop,ends.back()-starts[0])))/log(noStop));
      end = (beg + len) % G;
    }else if (index > 0){
      //Recombination interval begins at the start of other interval
      //Simulate the end-point via a truncated geometric distribution
      double r_3 = gsl_rng_uniform(rng);
      int len = floor(log(1-r_3*(1-pow(noStop,G+ends[index-1]-starts[index])))/log(noStop));
      end = (beg + len) % G;
    }else{
      //Recombination interval begins within an ancestral interval
      //Simulate the end-point via a truncated geometric distribution
      double r_3 = gsl_rng_uniform(rng);
      int len = floor(log(1-r_3*(1-pow(noStop,G-1)))/log(noStop));
      end = (beg + len) % G;
    }
  }
  //Check if end of recombinant interval falls between ancestral intervals
  if (end < starts[0]) end = ends.back();
  else{
    for (int a=1;a<b;++a){
      if ((end >= ends[a-1]) && (end < starts[a])){
        end == ends[a-1];
        break;
      }
    }
  }
}

void choose_clonalRecomb(const vector<double> &prob, const int G, const vector<int> &starts, const vector<int> &ends, int &beg, int &end, const double delta){
  //Choose a recombination interval for the given lineage which is clonal
  int b=starts.size();
  //Choose a start site of recombination at random
  double r_2=gsl_rng_uniform(rng);
  beg = 0;
  while (r_2 > prob[beg]){
    r_2 -= prob[beg];
    ++beg;
  }
  int len = gsl_ran_geometric(rng,1.0/delta);
  if (len > G-1){
    len = G-1;
  }
  end = (beg + len) % G;
  if ((starts[0] == 0) && (ends.back() == G-1)){
    if (b>1){
      //Check if end of recombinant interval falls between ancestral intervals
      for (int a=1;a<b;++a){
        if ((end >= ends[a-1]) && (end < starts[a])){
          end == ends[a-1];
          break;
        }
      }
    }
  }else{
    //Check if end of recombinant interval falls between ancestral intervals
    if (end < starts[0]) end = ends.back();
    else{
      for (int a=1;a<b;++a){
        if ((end >= ends[a-1]) && (end < starts[a])){
          end == ends[a-1];
          break;
        }
      }
    }
  }
}

void Arg::construct() {
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

  //Create ancestral material
  toCoalAncMat.push_back(vector<bool>(G,false));
  for (int ii=0;ii<b;++ii){
    for (int m=blockStarts[ii];m<blockEnds[ii];++m) toCoalAncMat.back()[m] = true;
  }

  //Calculate intervals for initial ancestral material
  intervalStarts.push_back(vector<int>(NULL));
  intervalEnds.push_back(vector<int>(NULL));
  calc_intervals(toCoalAncMat.back(), intervalStarts.back(), intervalEnds.back(), blockStarts, blockEnds, b);

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
      //Create initial ancestral material
      toCoalAncMat.push_back(toCoalAncMat[0]);

      //Set the probability of starting recombination at each site of the genome
      probStart.push_back(probStart[0]);
      recombRates.push_back(recombRates[0]);
      intervalStarts.push_back(intervalStarts[0]);
      intervalEnds.push_back(intervalEnds[0]);
  }

  t2=clock();
  cout << "Time spent intitiating nodes: " << t2-t1 << " which is: " << (double)(t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;

  double currentTime=0.0;
  double interval_time=0.0, recomb_probTime=0.0, recomb_intervalTime=0.0, combine_time=0.0, check_time=0.0;

  //Simulate the coalescence-recombination graph
  while (k>1) {
    //Calculate the current rate of recombination
    double currentRecomb = 0.0;
    for (int i=0;i<int(recombRates.size());++i) currentRecomb += recombRates[i];

    //Simulate time to next event exponentially
    currentTime+=gsl_ran_exponential(rng,2.0/(k*(k-1)+(2.0*currentRecomb)));

    //Randomly choose coalescence or recombination
    if (gsl_rng_uniform(rng)<1.0*(k*(k-1))/(k*(k-1)+(2.0*currentRecomb))){
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

      t1=clock();
      //Combine the ancestral material from the two children into the parent
      for (int m=0;m<b;++m){
        for (int a=blockStarts[m];a<blockEnds[m];++a){
          toCoalAncMat[i][a] = toCoalAncMat[i][a]||toCoalAncMat[j][a];
        }
      }
      t2=clock();
      combine_time += (t2-t1);
      t1=clock();
      //Test for fully coalesced material
      for (int m=0;m<b;++m){
        for (int a=blockStarts[m];a<blockEnds[m];++a) if (toCoalAncMat[i][a]==true){
          bool rem=true;
          for (int l=0;l<k;++l) if (l!=i&&toCoalAncMat[l][a]==true) {
            rem=false;
            break;
          }
          if (rem==true) toCoalAncMat[i][a]=false;
        }
      }
      t2=clock();
      check_time += t2-t1;
      intervalStarts[i] = vector<int>(NULL);
      intervalEnds[i] = vector<int>(NULL);

      t1 = clock();
      //Calculate the new sets of intervals of ancestral material for the coalesced node
      calc_intervals(toCoalAncMat[i], intervalStarts[i], intervalEnds[i], blockStarts, blockEnds, b);
      t2=clock();
      interval_time += t2-t1;
      t1 = clock();
      //Calculate the recombination rate for the new node
      if (clonal.back() == true){
        calc_clonalRecomb(G, delta, probStart[i], recombRates[i], intervalStarts[i], intervalEnds[i], noStop, siteRecomb);
      }else{
        calc_nonClonalRecomb(G, delta, probStart[i], recombRates[i], intervalStarts[i], intervalEnds[i], noStop, siteRecomb);
      }
      t2=clock();
      recomb_probTime += (t2-t1);
      //Remove the second child from the ARG
      toCoal[j] = toCoal.back();
      toCoal.pop_back();

      toCoalAncMat[j] = toCoalAncMat.back();
      toCoalAncMat.pop_back();

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

      //Recombinant interval does not wrap around the end of the genome
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
        if ((LTend >= blockStarts[m]) && (LTend <= blockStarts[m+1])){
          if (LTend >= blockEnds[m]) LTend = blocks[m+1]-1;
          else LTend = LTend - blockStarts[m] + blocks[m];
          break;
        }
      }
      changeLT[LTbeg] = true;
      if (LTend<(int)changeLT.size()) changeLT[LTend] = true;
      //Add the ages of the new nodes
      ages.push_back(currentTime);
      ages.push_back(currentTime);
      //First parent is clonal if the child is clonal
      clonal.push_back(clonal[toCoal[i]]);
      //Set second parent as non-clonal
      clonal.push_back(false);
      //Add the first node
      s.push_back(vector<int>(6,-1));
      //Add child of new parent
      s.back()[0]=toCoal[i];
      //Add second parent
      s.push_back(vector<int>(6,-1));
      //Add child of new parent
      s.back()[0]=toCoal[i];
      //Add start and end of import
      s.back()[4]=beg;
      s.back()[5]=end;
      //Set the parents of the child node
      s[toCoal[i]][2]=s.size()-2;
      s[toCoal[i]][3]=s.size()-1;
      //Set parent of child in the ARG
      toCoal[i]=s.size()-2;
      //Add non-clonal parent to the ARG
      toCoal.push_back(s.size()-1);
      //Set ancestral material of the new parents
      toCoalAncMat.push_back(vector<bool>(G,false));
      for (int ii=beg;ii<=end;++ii){
        toCoalAncMat.back()[ii]=toCoalAncMat[i][ii];
        toCoalAncMat[i][ii] = false;
      }
      
      toCoalAncMat.push_back(vector<bool>(G,false));
      for (int ii=blockStarts[0];ii <=end;++ii){
        toCoalAncMat.back()[ii]=toCoalAncMat[i][ii];
        toCoalAncMat[i][ii] = false;
      }
      for (int ii=beg;ii<blockEnds.back();++ii){
        toCoalAncMat.back()[ii]=toCoalAncMat[i][ii];
        toCoalAncMat[i][ii] = false;
      }

      intervalStarts[i] = vector<int>(NULL);
      intervalEnds[i] = vector<int>(NULL);
      t1=clock();
      //Calculate the intervals of ancestral material for the new nodes
      calc_intervals(toCoalAncMat[i], intervalStarts[i], intervalEnds[i], blockStarts, blockEnds, b);

      intervalStarts.push_back(vector<int>(NULL));
      intervalEnds.push_back(vector<int>(NULL));

      calc_intervals(toCoalAncMat.back(), intervalStarts.back(), intervalEnds.back(), blockStarts, blockEnds, b);
      t2=clock();
      interval_time += (t2-t1);

      //Calculate recombination rates and start-point probabilities for the new parents
      t1=clock();
      if (clonal[s.size()-2] == true){
        calc_clonalRecomb(G, delta, probStart[i], recombRates[i], intervalStarts[i], intervalEnds[i], noStop, siteRecomb);
      }else{
        calc_nonClonalRecomb(G, delta, probStart[i], recombRates[i], intervalStarts[i], intervalEnds[i], noStop, siteRecomb);
      }

      //And for non-clonal parent
      recombRates.push_back(0.0);
      probStart.push_back(vector<double>(G,0.0));
      if (clonal.back() == true){
        calc_clonalRecomb(G, delta, probStart.back(), recombRates.back(), intervalStarts.back(), intervalEnds.back(), noStop, siteRecomb);
      }else{
        calc_nonClonalRecomb(G, delta, probStart.back(), recombRates.back(), intervalStarts.back(), intervalEnds.back(), noStop, siteRecomb);
      }
      t2=clock();
      recomb_probTime += (t2-t1);

      cout << "Recomb interval: " << beg << " " << end << endl;
      for (int a=0;a<intervalStarts[i].size();++a) cout << intervalStarts[i][a] << " " << intervalEnds[i][a] << " ";
      cout << endl;
      for (int a=0;a<intervalStarts.back().size();++a) cout << intervalStarts.back()[a] << " " << intervalEnds.back()[a] << " ";
      cout << endl;
      cout << endl;

      ++k;
    }
  }
  cout << "Time spent on combining ancestries: " << combine_time << " which is : " << (double)combine_time/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent on calculating intervals: " << interval_time << " which is : " << (double)interval_time/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent on finding recombinant intervals: " << recomb_intervalTime << " which is : " << (double)recomb_intervalTime/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent on calculating recomb probabilities: " << recomb_probTime << " which is : " << (double)recomb_probTime/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time spent checking for coalescence: " << check_time << " which is : " << (double)check_time/CLOCKS_PER_SEC << " seconds" << endl;
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
          if (beg <= nd){
            for (int j=beg;j<=nd;j++) genotypes[i]->at(j)=genotypes[s[i][3]]->at(j);
          }else{
            //TB - interval wraps around genome
            for (int j=0;j<=nd;++j) genotypes[i]->at(j)=genotypes[s[i][3]]->at(j);
            for (int j=beg;j<L;++j) genotypes[i]->at(j)=genotypes[s[i][3]]->at(j);
          }
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
        if (beg <= nd){
          if (site>=beg && site<=nd) s4[s[k][3]]=true;else s4[s[k][2]]=true;
        }else{
          //TB - interval wraps around genome
          if (site<=nd || site>=beg) s4[s[k][3]]=true;else s4[s[k][2]]=true;
        }
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
    } else {
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
        } else if (s[s[i][0]][2]==i){ //recipient
          if (s[s[s[i][0]][3]][4] <= s[s[s[i][0]][3]][5]){
            for (int k=0;k<L;++k) {if (k>=s[s[s[i][0]][3]][4]&&k<=s[s[s[i][0]][3]][5]) ancmat.back()[k]=false;}
          }else{
            //TB - recombinant interval wraps around the genome
            for (int k=0;k<L;++k){
              if(k>=s[s[s[i][0]][3]][4])ancmat.back()[k]=false;
              else if (k<=s[s[s[i][0]][3]][5]) ancmat.back()[k]=false;
            }
          }
        }else{//donor
          if (s[i][4] <= s[i][5]){
            for (int k=0;k<L;++k) {if (k<s[i][4]||k>s[i][5]) ancmat.back()[k]=false;}
          }else{
            for (int k=0;k<L;++k) {if (k<s[i][5]&&k<s[i][4]) ancmat.back()[k]=false;}
          }
        }
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
