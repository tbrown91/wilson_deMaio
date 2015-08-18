#ifndef MRCA_CHECKS
#define MRCA_CHECKS

void removeAncMat(const int start, const int end, list<int> &starts, list<int> &ends){
  //Reove MRCA material from chosen nodes
  list<int>::iterator itStart=starts.begin(), itEnd = ends.begin();
  while (itEnd!=ends.end()){
    if (end < *itStart) return; //No more ancestral material to check
    else if (start <= *itStart){
      if (end >= *itEnd){
        //Entire interval to be removed
        itStart = starts.erase(itStart);
        itEnd = ends.erase(itEnd);
      }else{
        //Start of interval to be removed
        *itStart = end + 1;
        return; //No more ancestral material to check
      }
    }else if (start <= *itEnd){
      if (end >= *itEnd){
        //End of interval to be removed
        *itEnd = start - 1;
        ++itStart;
        ++itEnd;
      }else{
        //Middle of interval to be removed
        ends.insert(itEnd,start - 1);
        ++itStart;
        starts.insert(itStart,end + 1);
        return;//No more material to check;
      }
    }else{
        ++itStart;
        ++itEnd;
    }
  }
}

void modifyMRCA(Arg::MRCA &M, const int start, const int end){
//update the list of MRCA intervals by subtracting an interval overlapping in the two coalesceing lineages
  M.itStart=M.starts.begin();
  M.itEnd=M.ends.begin();
  M.itValue=M.values.begin();
	if (M.itStart==M.starts.end()){
		cout << "Error : this part should be in the MRCA vector, but it isn't (1) " << *M.itStart << " " << *M.starts.end() << endl;
	 	return;
	}
	while (*(M.itEnd) < start) {
		++(M.itStart);
		++(M.itEnd);
		++(M.itValue);
		if (M.itStart==M.starts.end()){
			cout << "Error : this part should be in the MRCA vector, but it isn't (2)" << *M.itStart << " " << *M.starts.end() << endl;
      break;
		}
	}
	if (start<*(M.itStart)){
		cout << "Error : this part should be in the MRCA vector, but it isn't (3)" << start << " " << *M.itStart << endl;
    int a;
    cin >> a;
    return;
	}
	//first element in MRCA list that overlaps, but starting values do not coincide
	if ((M.itStart!=M.starts.end()) && (*M.itStart<start)){
    ++(M.itStart);
    M.starts.insert(M.itStart,start);
    --(M.itStart);
		M.ends.insert(M.itEnd,start-1);
		M.values.insert(M.itValue,*M.itValue);
	}
	//iteratively look at all intervals in the MRCA structure overlapping completely the given interval
	while ((M.itStart!=M.starts.end()) && (*M.itEnd<=end)) {
    --*(M.itValue);
		++(M.itStart);
		++(M.itEnd);
		++(M.itValue);
	}
  // look at the final interval in the MRCA, if present, that overlaps the given interval but has end after the end of the given interval
  if ((M.itStart!=M.starts.end()) && (*M.itStart<=end)){
    ++(M.itStart);
    M.starts.insert(M.itStart,end+1);
    M.ends.insert(M.itEnd,end);
    M.values.insert(M.itValue,*M.itValue);
    --M.itValue;
    --*(M.itValue);
  }
}

void update_MRCA(Arg::MRCA &M, list<int> &starts_1, list<int> &ends_1, const list<int> &starts_2, const list<int> &ends_2){
  //Find overlapping regions of ancestral material in the two children and update the MRCA lists
  if ((starts_1.size() == 0) || (starts_2.size() == 0)) return; //No ancestral material in one or both nodes, therefore MRCA struct does not need to be updated
  list<int>::iterator itStart1 = starts_1.begin(), itEnd1 = ends_1.begin();
  list<int>::const_iterator itStart2 = starts_2.begin(), itEnd2 = ends_2.begin();
  int currentStart1=0, currentStart2=0;
  currentStart1=*itStart1, currentStart2=*itStart2;
  while ((itStart1 != starts_1.end()) && (itStart2 != starts_2.end())){
    if ((currentStart1 < currentStart2) && (*itEnd1>=currentStart2)){ //overlap
      currentStart1=currentStart2;
    }else if ((currentStart2 < currentStart1) && (*itEnd2>=currentStart1)){ //overlap
      currentStart2=currentStart1;
    }else if ((currentStart1 < currentStart2) && (*itEnd1<currentStart2)){ //no overlap
      ++itStart1;
      ++itEnd1;
      if (itStart1!=starts_1.end()) currentStart1=*itStart1;
    }else if ((currentStart2 < currentStart1) && (*itEnd2<currentStart1)){ //no overlap
      ++itStart2;
      ++itEnd2;
      if (itStart2 != starts_2.end()) currentStart2=*itStart2;
    }else if ((currentStart1 == currentStart2) && (*itEnd1<*itEnd2)){ //overlap, same start
      modifyMRCA(M, currentStart1, *itEnd1);
      currentStart2=*itEnd1+1;
      ++itStart1;
      ++itEnd1;
      if (itStart1 != starts_1.end()) currentStart1=*itStart1;
    }else if ((currentStart1 == currentStart2) && (*itEnd2<*itEnd1)){ //overlap, same start
      modifyMRCA(M, currentStart1, *itEnd2);
      currentStart1=*itEnd2+1;
      ++itStart2;
      ++itEnd2;
      if (itStart2 != starts_2.end()) currentStart2=*itStart2;
    }else if ((currentStart1 == currentStart2) && (*itEnd2==*itEnd1)){ //overlap, same start, same end
      modifyMRCA(M, currentStart1, *itEnd1);
      ++itStart2;
      ++itEnd2;
      if (itStart2 != starts_2.end()) currentStart2=*itStart2;
      ++itStart1;
      ++itEnd1;
      if (itStart1 != starts_1.end()) currentStart1=*itStart1;
    }else{
      cout << "Error : this case should not happen in updating MRCA vector" << endl;
      return;
    }
  }
}
#endif
