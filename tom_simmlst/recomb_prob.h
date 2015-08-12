#include <math.h>

void calc_nonClonalRecomb(const int G, const double delta, vector<double> &prob, double &recombRate, const vector<int> &starts, const vector<int> &ends, const double noStop, const double siteRecomb, int &totMaterial){
  //Total amount of ancestral material in the node
  totMaterial = 0;
  int b = starts.size();
  //Reset recombination rate and probability of recombination along the genome
  prob = vector<double>(b,0.0);
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
        prob[a] = siteRecomb*delta*(1-pow(noStop,starts[a]-ends[a-1]))*(1-pow(noStop,G-starts[a]+ends[a-1]));
        recombRate += prob[a];
      }
      recombRate += siteRecomb*(totMaterial-(b-1))*(1-pow(noStop,G-1));
      prob[0] = siteRecomb*(1-pow(noStop,G-1));
    }else{
      //Interval is the entire genome
      recombRate = siteRecomb*totMaterial*(1-pow(noStop,G-1));
      prob[0] = siteRecomb*(1-pow(noStop,G-1));
    }
  }else{
    //Intervals are contained within the genome
    prob[0] = siteRecomb*delta*(1-pow(noStop,G+starts[0]-ends.back()))*(1-pow(noStop,ends.back()-starts[0]));
    recombRate += prob[0];
    for (int a=1;a<b;++a){
      prob[a] = siteRecomb*delta*(1-pow(noStop,starts[a]-ends[a-1]))*(1-pow(noStop,G-starts[a]+ends[a-1]));
      recombRate += prob[a];
    }
    recombRate += siteRecomb*(totMaterial-b)*(1-pow(noStop,G-1));
  }
}

void calc_clonalRecomb(const int G, const double delta, vector<double> &prob, double &recombRate, const vector<int> &starts, const vector<int> &ends, const double noStop, const double siteRecomb, int &totMaterial){
  //Total amount of ancestral material in the node
  totMaterial = 0;
  int b = starts.size();
  //Reset recombination rate and probability of recombination along the genome
  prob = vector<double>(b,0.0);
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
        prob[a] = siteRecomb*delta*(1-pow(noStop,starts[a]-ends[a-1]));
        recombRate += prob[a];
      }
      recombRate += siteRecomb*(totMaterial-(b-1));
      prob[0] = siteRecomb;
    }else{
      //Interval is the entire genome
      recombRate = siteRecomb*totMaterial;
      prob[0] = siteRecomb;
    }
  }else{
    //Intervals are contained within the genome
    prob[0] = siteRecomb*delta*(1-pow(noStop,G+starts[0]-ends.back()));
    recombRate += prob[0];
    for (int a=1;a<b;++a){
      prob[a] = siteRecomb*delta*(1-pow(noStop,starts[a]-ends[a-1]));
      recombRate += prob[a];
    }
    recombRate += siteRecomb*(totMaterial-b);
  }
}
