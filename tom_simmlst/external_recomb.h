#ifndef EXTERNAL_RECOMB
#define EXTERNAL_RECOMB

void external_interval(const int beg, const int end, int &recombStart, int &recombEnd, const list<int> &starts, const list<int> &ends, int &ext_check){
  //Choose an interval for the external recombinant event
  list<int>::const_iterator itS = starts.begin(), itE = ends.begin();
  if (beg <= end){
    //Recombinant interval does not wrap around end of genome
    if ((end < starts.front()) || (beg > ends.back())){//Interval falls outside of ancestral material
      ext_check = 1;
      return;
    }
    ++itS;
    while (itS != starts.end()){
      if ((beg > *itE) && (end < *itS)){//Interval falls outside of ancestral material
        ext_check = 1;
        return;
      }
      ++itS;
      ++itE;
    }
    itS = starts.begin(), itE = ends.begin();
    while (itS != starts.end()){
      if (beg <= *itS) recombStart = *itS;
      else if (beg <= *itE) recombStart = beg;
      if (end >= *itE) recombEnd = *itE;
      else if (end >= *itS) recombEnd = end;
      else break;
      ++itS;
      ++itE;
    }
  }else{
    if ((beg > ends.back()) && (end < starts.front())){//Interval falls outside of ancestral material
      ext_check = 1;
      return;
    }
    ++itS;
    while (itS != starts.end()){
      if ((end > *itE) && (beg < *itS)){//Interval falls outside of ancestral material
        ext_check = 1;
        return;
      }
      ++itS;
      ++itE;
    }
    itS = starts.begin(), itE = starts.begin();
    if (end <= *itS) recombEnd = ends.back();
    if (beg >= ends.back()) recombStart = ends.back();
    ++itS;
    while (itS != starts.end()){
      if (end <= *itE){
        recombEnd = end;
        break;
      }else if (end < *itS){
        recombEnd = *itE;
        break;
      }
      ++itS;
      ++itE;
    }
    itS = starts.begin(), itE = ends.begin();
    while (itS != starts.end()){
      if (beg <= *itS){
        recombStart = *itS;
        break;
      }else if (beg <= *itE){
        recombStart = beg;
        break;
      }
      ++itS;
      ++itE;
    }
  }
}
#endif
