#include <math.h>

void calc_nonClonalRecomb(const int G, const double delta, vector<double> &prob, double &recombRate, const vector<int> &starts, const vector<int> &ends, const double noStop, const double siteRecomb){
  //Total amount of ancestral material in the node
  int totMaterial = 0;
  int b = starts.size();
  //Reset recombination rate and probability of recombination along the genome
  prob = vector<double>(G,0.0);
  recombRate = 0;
  if (b == 0) return;
  for (int m=0;m<b;++m){
    totMaterial += (ends[m] - starts[m]) + 1;
  }
  if (totMaterial <= 1) return;
  if ((starts[0] == 0) && (ends.back() == G-1)){
    //Interval wraps around the end of the genome
    if (b > 1){
      for (int a=1;a<b;++a){
        prob[starts[a]] = siteRecomb*delta*(1-pow(noStop,starts[a]-ends[a-1]))*(1-pow(noStop,G-starts[a]+ends[a-1]));
        recombRate += prob[starts[a]];
      }
      recombRate += siteRecomb*(totMaterial-(b-1))*(1-pow(noStop,G-1));
      double recombProb = siteRecomb*(1-pow(noStop,G-1))/recombRate;
      for (int a=starts[0];a<=ends[0];++a) prob[a] = recombProb;
      for (int ii=1;ii<b;++ii){
        for (int a=starts[ii]+1;a<=ends[ii];++a) prob[a] = recombProb;
        prob[starts[ii]] /= recombRate;
      }
    }else{
      //Interval is the entire genome
      recombRate = siteRecomb*totMaterial*(1-pow(noStop,G-1));
      double recombProb = 1/(double)totMaterial;
      for (int a=0;a<G;++a) prob[a] = recombProb;
    }
  }else{
    //Intervals are contained within the genome
    prob[starts[0]] = siteRecomb*delta*(1-pow(noStop,G+starts[0]-ends.back()))*(1-pow(noStop,ends.back()-starts[0]));
    recombRate += prob[starts[0]];
    for (int a=1;a<b;++a){
      prob[starts[a]] = siteRecomb*delta*(1-pow(noStop,starts[a]-ends[a-1]))*(1-pow(noStop,G-starts[a]+ends[a-1]));
      recombRate += prob[starts[a]];
    }
    recombRate += siteRecomb*(totMaterial-b)*(1-pow(noStop,G-1));
    double recombProb = siteRecomb*(1-pow(noStop,G-1))/recombRate;
    for (int ii=0;ii<b;++ii){
      for (int a=starts[ii]+1;a<=ends[ii];++a) prob[a] = recombProb;
      prob[starts[ii]] /= recombRate;
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
  if (b == 0) return;
  for (int m=0;m<b;++m){
    totMaterial += (ends[m] - starts[m]) + 1;
  }
  if (totMaterial <= 1) return;
  if ((starts[0] == 0) && (ends.back() == G-1)){
    //Interval wraps around the end of the genome
    if (b > 1){
      for (int a=1;a<b;++a){
        prob[starts[a]] = siteRecomb*delta*(1-pow(noStop,starts[a]-ends[a-1]));
        recombRate += prob[starts[a]];
      }
      recombRate += siteRecomb*(totMaterial-(b-1));
      double recombProb = siteRecomb/recombRate;
      for (int a=0;a<=ends[0];++a) prob[a] = recombProb;
      for (int ii=1;ii<b;++ii){
        for (int a=starts[ii]+1;a<=ends[ii];++a) prob[a] = recombProb;
        prob[starts[ii]] /= recombRate;
      }
    }else{
      //Interval is the entire genome
      recombRate = siteRecomb*totMaterial;
      double recombProb = 1/(double)totMaterial;
      for (int a=0;a<G;++a) prob[a] = recombProb;
    }
  }else{
    //Intervals are contained within the genome
    prob[starts[0]] = siteRecomb*delta*(1-pow(noStop,G+starts[0]-ends.back()));
    recombRate += prob[starts[0]];
    for (int a=1;a<b;++a){
      prob[starts[a]] = siteRecomb*delta*(1-pow(noStop,starts[a]-ends[a-1]));
      recombRate += prob[starts[a]];
    }
    recombRate += siteRecomb*(totMaterial-b);
    double recombProb = siteRecomb/recombRate;
    for (int ii=0;ii<b;++ii){
      for (int a=starts[ii]+1;a<=ends[ii];++a) prob[a] = recombProb;
      prob[starts[ii]] /= recombRate;
    }
  }
}
