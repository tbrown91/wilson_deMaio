//update the list of MRCA intervals by subtracting an interval overlapping in the two coalesceing lineages
void modifyMRCA(MRCA &M, int start, int end){
	if ((*M).itStart==(*M).starts.end()){
		cout << "Error : this part should be in the MRCA vector, but it isn't (1)" << endl;
    	exit();
	}
	while (*((*M).itEnd) < start) {
		++((*M).itStart);
		++((*M).itEnd);
		++((*M).itValue);
		if ((*M).itStart==((*M).starts).end()){
			cout << "Error : this part should be in the MRCA vector, but it isn't (2)" << endl;
    		//exit();
      break;
		}
	}

	if (start<*((*M).itStart)){
		cout << "Error : this part should be in the MRCA vector, but it isn't (3)" << endl;
    	//exit();
    break;
	}

	//first element in MRCA list that overlaps, but starting values do not coincide
	if (*((*M).itStart)<start){
		(*M).starts.insert((*M).itStart,start);
		(*M).ends.insert((*M).itEnd,start-1);
		(*M).values.insert((*M).itValue,*((*M).itValue));
		++((*M).itStart);
		++((*M).itEnd);
		++((*M).itValue);
		*((*M).itStart)=start;
	}

	//iteratively look at all intervals in the MRCA structure overlapping completely the given interval
	while (((*M).itStart!=(*M).starts.end()) && (*((*M).itEnd)<=end)) {
		*((*M).itValue)=*((*M).itValue)-1;
		if (*((*M).itValue)==0) {
      cout << "Error : MRCA values should not reach 0!" << endl;
    		//exit();
      return;
		}else if (*((*M).itValue)==1) {
			//removeFromAncestralMaterial(*((*M).itStart),*((*M).itEnd));/////////////////////////////////An interval has just reached MRCA. What do you want to do with this information?

		}
		++((*M).itStart);
		++((*M).itEnd);
		++((*M).itValue);
	}

	// look at the final interval in the MRCA, if present, that overlaps the given interval but has end after the end of the given interval
	if (((*M).itStart!=(*M).starts.end()) && (*((*M).itStart)<=end)){
		(*M).starts.insert((*M).itStart,*((*M).itStart));
		(*M).ends.insert((*M).itEnd,end);
		(*M).values.insert((*M).itValue,*((*M).itValue)-1);
		++((*M).itStart);
		++((*M).itEnd);
		++((*M).itValue);
		*((*M).itStart)=end+1;
	}
}

//Update MRCA structure providing intervals that are shared by the two coalescing lineages
void update_MRCA(MRCA &M, vector<int> &starts1, vector<int> &ends1, vector<int> &starts2, vector<int> &ends2){
  int index1=0, index2=0;
  int b1 = starts1.size(), b2 = starts2.size();
  (*M).itStart=(*M).starts.begin();
  (*M).itEnd=(*M).ends.begin();
  (*M).itValue=(*M).values.begin();
  int currentStart1=0, currentStart2=0;
  if ((index1 != b1) && (index2 != b2)) currentStart1=start1[index1], currentStart2=start2[index2];
  while ((index1 != b1) && (index2 != b2)){
    if ((currentStart1 < currentStart2) && (ends1[index1]>=currentStart2)){ //overlap
    	currentStart1=currentStart2;
    }else if ((currentStart2 < currentStart1) && (ends2[index2]>=currentStart1)){ //overlap
    	currentStart2=currentStart1;
    }else if ((currentStart1 < currentStart2) && (ends1[index1]<currentStart2)){ //no overlap
    	++index1;
    	if (index1<b1) currentStart1=starts1[index1];
    }else if ((currentStart2 < currentStart1) && (ends2[index2]<currentStart1)){ //no overlap
    	++index2;
    	if (index2<b2) currentStart2=starts2[index2];
    }else if ((currentStart1 == currentStart2) && (ends1[index1]<ends2[index2])){ //overlap, same start
    	modifyMRCA(M, currentStart1, ends1[index1]);
    	currentStart2=ends1[index1]+1;
    	++index1;
    	if (index1<b1) currentStart1=starts1[index1];
    }else if ((currentStart1 == currentStart2) && (ends2[index2]<ends1[index1])){ //overlap, same start
    	modifyMRCA(M, currentStart2, ends2[index2]);
    	currentStart1=ends2[index2]+1;
    	++index2;
    	if (index2<b2) currentStart2=starts2[index2];
    }else if ((currentStart1 == currentStart2) && (ends2[index2]==ends1[index1])){ //overlap, same start, same end
    	modifyMRCA(M, currentStart2, ends2[index2]);
    	++index2;
    	if (index2<b2) currentStart2=starts2[index2];
    	++index1;
    	if (index1<b1) currentStart1=starts1[index1];
    }else{
    	cout << "Error : this case should not happen in updating MRCA vector" << endl;
    	//exit();
      break;
    }
  }
}
