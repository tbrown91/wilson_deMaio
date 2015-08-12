#ifndef MODIFY_MRCA
#define MODIFY_MRCA
void modifyMRCA(MRCA& M, int start, int end, vector<int> &starts1, vector<int> &ends1, vector<int> &starts2, vector<int> &ends2){
//update the list of MRCA intervals by subtracting an interval overlapping in the two coalesceing lineages
	if ((M).itStart==(M).starts.end()){
		cout << "Error : this part should be in the MRCA vector, but it isn't (1)" << endl;
		return;
    //exit();
	}
	while (*((M).itEnd) < start) {
		++((M).itStart);
		++((M).itEnd);
		++((M).itValue);
		if ((M).itStart==((M).starts).end()){
			cout << "Error : this part should be in the MRCA vector, but it isn't (2)" << endl;
    	//exit();
      break;
		}
	}

	if (start<*((M).itStart)){
		cout << "Error : this part should be in the MRCA vector, but it isn't (3)" << endl;
  	//exit();
    break;
	}

	//first element in MRCA list that overlaps, but starting values do not coincide
	if (*((M).itStart)<start){
		(M).starts.insert((M).itStart,start);
		(M).ends.insert((M).itEnd,start-1);
		(M).values.insert((M).itValue,*((M).itValue));
		++((M).itStart);
		++((M).itEnd);
		++((M).itValue);
		*((M).itStart)=start;
	}

	//iteratively look at all intervals in the MRCA structure overlapping completely the given interval
	while (((M).itStart!=(M).starts.end()) && (*((M).itEnd)<=end)) {
		*((M).itValue)=*((M).itValue)-1;
		if (*((M).itValue)==0) {
      cout << "Error : MRCA values should not reach 0!" << endl;
    	//exit();
      return;
		}else if (*((M).itValue)==1) {
			//removeFromAncestralMaterial(*((M).itStart),*((M).itEnd));/////////////////////////////////An interval has just reached MRCA. What do you want to do with this information?
			return;
		}
		++((M).itStart);
		++((M).itEnd);
		++((M).itValue);
	}

	// look at the final interval in the MRCA, if present, that overlaps the given interval but has end after the end of the given interval
	if (((M).itStart!=(M).starts.end()) && (*((M).itStart)<=end)){
		(M).starts.insert((M).itStart,*((M).itStart));
		(M).ends.insert((M).itEnd,end);
		(M).values.insert((M).itValue,*((M).itValue)-1);
		++((M).itStart);
		++((M).itEnd);
		++((M).itValue);
		*((M).itStart)=end+1;
	}
}
#endif
