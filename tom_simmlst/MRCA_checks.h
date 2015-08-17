#ifndef MRCA_CHECKS
#define MRCA_CHECKS

void removeAncMat(const int start, const int end, vector<int> &starts, vector<int> &ends){
  //Reove MRCA material from chosen nodes
  vector<int> tempStarts;
  vector<int> tempEnds;
  for (size_t a=0;a<starts.size();++a){
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

void update_MRCA(Arg::MRCA &M, vector<int> &starts_1, vector<int> &ends_1, vector<int> &starts_2, vector<int> &ends_2){
  //Find overlapping regions of ancestral material in the two children and update the MRCA lists
  int index1=0, index2=0;
  int b1 = starts_1.size(), b2 = starts_2.size();
  int currentStart1=0, currentStart2=0;
  if ((index1 != b1) && (index2 != b2)) currentStart1=starts_1[index1], currentStart2=starts_2[index2];
  while ((index1 != b1) && (index2 != b2)){
    if ((currentStart1 < currentStart2) && (ends_1[index1]>=currentStart2)){ //overlap
      currentStart1=currentStart2;
    }else if ((currentStart2 < currentStart1) && (ends_2[index2]>=currentStart1)){ //overlap
      currentStart2=currentStart1;
    }else if ((currentStart1 < currentStart2) && (ends_1[index1]<currentStart2)){ //no overlap
      ++index1;
      if (index1<b1) currentStart1=starts_1[index1];
    }else if ((currentStart2 < currentStart1) && (ends_2[index2]<currentStart1)){ //no overlap
      ++index2;
      if (index2<b2) currentStart2=starts_2[index2];
    }else if ((currentStart1 == currentStart2) && (ends_1[index1]<ends_2[index2])){ //overlap, same start
      modifyMRCA(M, currentStart1, ends_1[index1]);
      currentStart2=ends_1[index1]+1;
      ++index1;
      if (index1<b1) currentStart1=starts_1[index1];
    }else if ((currentStart1 == currentStart2) && (ends_2[index2]<ends_1[index1])){ //overlap, same start
      modifyMRCA(M, currentStart1, ends_2[index2]);
      currentStart1=ends_2[index2]+1;
      ++index2;
      if (index2<b2) currentStart2=starts_2[index2];
    }else if ((currentStart1 == currentStart2) && (ends_2[index2]==ends_1[index1])){ //overlap, same start, same end
      modifyMRCA(M, currentStart1, ends_1[index1]);
      ++index2;
      if (index2<b2) currentStart2=starts_2[index2];
      ++index1;
      if (index1<b1) currentStart1=starts_1[index1];
    }else{
      cout << "Error : this case should not happen in updating MRCA vector" << endl;
      break;
    }
  }
}
#endif
