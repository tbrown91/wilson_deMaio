#ifndef RECOMB_EVENT
#define RECOMB_EVENT
#include <math.h>

void split_ancestries(vector<int> &starts_1, vector<int> &ends_1, vector<int> &starts_2, vector<int> &ends_2, const int beg, const int end){
  //Split the ancestry of the recombinant node into two, extracting the recombinant interval
  int b = starts_1.size();
  vector<int> tempStarts_1;
  vector<int> tempEnds_1;
  if (beg <= end){
    //Recombinant break does not wrap around the end of the genome
    for (int a=0;a<b;++a){
      //Find the interval containing the start and end points of the recombinant interval
      if ((beg > ends_1[a]) || (end < starts_1[a])){
        //Ancestral interval is not in recombinant break
        tempStarts_1.push_back(starts_1[a]);
        tempEnds_1.push_back(ends_1[a]);
      }else if((beg <= starts_1[a]) && (end >= ends_1[a])){
        //Recombinant break takes entire ancestral interval
        starts_2.push_back(starts_1[a]);
        ends_2.push_back(ends_1[a]);
      }else if((beg <= starts_1[a]) && (end >= starts_1[a])){
        //Break falls at start of ancestral interval
        starts_2.push_back(starts_1[a]);
        ends_2.push_back(end);
        tempStarts_1.push_back(end+1);
        tempEnds_1.push_back(ends_1[a]);
      }else if((beg <= ends_1[a]) && (end >= ends_1[a])){
        //Break falls at end of ancestral interval
        tempStarts_1.push_back(starts_1[a]);
        tempEnds_1.push_back(beg-1);
        starts_2.push_back(beg);
        ends_2.push_back(ends_1[a]);
      }else{
        //Break falls inside ancestral interval
        tempStarts_1.push_back(starts_1[a]);
        tempEnds_1.push_back(beg-1);
        starts_2.push_back(beg);
        ends_2.push_back(end);
        tempStarts_1.push_back(end+1);
        tempEnds_1.push_back(ends_1[a]);
      }
    }
  }else{
    //Recombinant break wraps around the end of the genome
    if (beg == (end+1)){
      //Recombinant interval takes entire genome
      starts_2 = starts_1;
      ends_2 = ends_1;
    }
    else{
      for (int a=0;a<b;++a){
        if (end >= ends_1[a]){
          //Break covers the ancestral inerval at the start
          starts_2.push_back(starts_1[a]);
          ends_2.push_back(ends_1[a]);
        }else if (beg <= starts_1[a]){
          //Break covers the ancestral inerval at the end
          starts_2.push_back(starts_1[a]);
          ends_2.push_back(ends_1[a]);
        }else if ((end >= starts_1[a]) && (beg <= ends_1[a])){
          //Break covers start and end of interval
          starts_2.push_back(starts_1[a]);
          ends_2.push_back(end);
          tempStarts_1.push_back(end+1);
          tempEnds_1.push_back(beg-1);
          starts_2.push_back(beg);
          ends_2.push_back(ends_1[a]);
        }else if (((end >= starts_1[a]) && (end < ends_1[a])) && (beg > ends_1[a])){
          //Break falls across an interval at the start
          starts_2.push_back(starts_1[a]);
          ends_2.push_back(end);
          tempStarts_1.push_back(end+1);
          tempEnds_1.push_back(ends_1[a]);
        }else if (((beg <= ends_1[a]) && (beg > starts_1[a])) && (end < starts_1[a])){
          //Break falls across an interval at the end
          tempStarts_1.push_back(starts_1[a]);
          tempEnds_1.push_back(beg-1);
          starts_2.push_back(beg);
          ends_2.push_back(ends_1[a]);
        }else{
          //Break does not cover ancestral inerval
          tempStarts_1.push_back(starts_1[a]);
          tempEnds_1.push_back(ends_1[a]);
        }
      }
    }
  }
  starts_1 = tempStarts_1;
  ends_1 = tempEnds_1;
}

void choose_nonClonalRecomb(const vector<double> &prob, const int G, const vector<int> &starts, const vector<int> &ends, int &beg, int &end, const double noStop, const int totMaterial, const double recombRate){
  //Choose a recombination interval for the chosen lineage which is non-clonal
  int b=starts.size();
  double r_1 = (gsl_rng_uniform(rng)*recombRate);
  int index = 0;
  while (r_1 > prob[index]){
    r_1 -= prob[index];
    ++index;
    if (index == b) break;
  }
  if (index == b){
    //Choose a start point within the ancestral intervals
    beg = (int)floor(gsl_rng_uniform(rng)*(totMaterial-b));
    beg += starts[0]+1;
    for (int a=0;a<b;++a){
      if (beg > ends[a]){
        beg = beg - ends[a] + starts[a+1];
      }else break;
    }
    double r_2 = gsl_rng_uniform(rng);
    int len = (int)floor(log(1-r_2*(1-pow(noStop,G-1)))/log(noStop));
    end = (beg + len) % G;
  }else{
    //Start of recombination occurs at beginning of an ancestral interval
    if (index == 0){
      if ((starts[0] == 0) && (ends.back() == G-1)){
        //Start site is at beginning of genome inside an ancestral interval
        beg = 0;
        double r_2 = gsl_rng_uniform(rng);
        int len = (int)floor(log(1-r_2*(1-pow(noStop,G-1)))/log(noStop));
        end = (beg + len) % G;
      }else{
        //Interval starts at first interval
        beg = starts[0];
        double r_2 = gsl_rng_uniform(rng);
        int len = (int)floor(log(1-r_2*(1-pow(noStop,ends.back()-starts[0])))/log(noStop));
        end = (beg + len) % G;
      }
    }else{
      beg = starts[index];
      //Simulate recombinant break length via a truncated geometric distribution
      double r_2 = gsl_rng_uniform(rng);
      int len = (int)floor(log(1-r_2*(1-pow(noStop,G+ends[index-1]-starts[index])))/log(noStop));
      end = (beg + len) % G;
    }
  }
  if ((starts[0] == 0) && (ends.back() == G-1)){
    //Ancestral material wraps around end of genome
    if (b>1){
      //Check if end of recombinant interval falls between ancestral intervals
      for (int a=1;a<b;++a){
        if ((end >= ends[a-1]) && (end < starts[a])){
          end = ends[a-1];
          break;
        }else if (end < ends[a-1]) break;
      }
    }
  }else{
    //Check if end of recombinant interval falls between ancestral intervals
    if ((end < starts[0]) || (end > ends.back())) end = ends.back();
    else{
      for (int a=1;a<b;++a){
        if ((end >= ends[a-1]) && (end < starts[a])){
          end = ends[a-1];
          break;
        }else if (end < ends[a-1]) break;
      }
    }
  }
}

void choose_clonalRecomb(const vector<double> &prob, const int G, const vector<int> &starts, const vector<int> &ends, int &beg, int &end, const double delta, const int totMaterial, const double recombRate){
  //Choose a recombination interval for the chosen lineage which is non-clonal
  int b=starts.size();
  double r_1 = (gsl_rng_uniform(rng)*recombRate);
  int index = 0;
  while (r_1 > prob[index]){
    r_1 -= prob[index];
    ++index;
    if (index == b) break;
  }
  if (index == b){
    //Choose a start point within the ancestral intervals
    beg = (int)floor(gsl_rng_uniform(rng)*(totMaterial-b));
    beg += starts[0]+1;
    for (int a=0;a<b;++a){
      if (beg > ends[a]){
        beg = beg - ends[a] + starts[a+1];
      }else break;
    }
  }else{
    beg = starts[index];
  }
  int len = gsl_ran_geometric(rng,1.0/delta);
  if (len > G-1){
    len = G-1;
  }
  end = (beg + len) % G;
  if ((starts[0] == 0) && (ends.back() == G-1)){
    //Ancestral material wraps around end of genome
    if (b>1){
      //Check if end of recombinant interval falls between ancestral intervals
      for (int a=1;a<b;++a){
        if ((end >= ends[a-1]) && (end < starts[a])){
          end = ends[a-1];
          break;
        }else if (end < ends[a-1]) break;
      }
    }
  }else{
    //Check if end of recombinant interval falls between ancestral intervals
    if ((end < starts[0]) || (end > ends.back())) end = ends.back();
    else{
      for (int a=1;a<b;++a){
        if ((end >= ends[a-1]) && (end < starts[a])){
          end = ends[a-1];
          break;
        }else if (end < ends[a-1]) break;
      }
    }
  }
}
#endif
