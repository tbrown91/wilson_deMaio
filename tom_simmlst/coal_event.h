#ifndef COAL_EVENT
#define COAL_EVENT
// add an interval to a list, or join it to the last interval of the list if they overlap
void add_interval_to_interval_list(list<int> &starts, list<int> &ends, int start, int end){
	list<int>::iterator itEnd = ends.end();
	--itEnd;
	if (starts.size() == 0){
		//Add the first interval to the list
		starts.push_back(start);
    ends.push_back(end);
  }else{
		//If the new interal overlaps with the old interval, increase the length of the old interval
  	if (start<=(*itEnd+1)){
  		if (end>*itEnd){
  			*itEnd=end;
  		}
  	}else{
			//Add the new interval to the list
  		starts.push_back(start);
  		ends.push_back(end);
  	}
  }
}

void combine_ancestries(list<int> &starts_1, list<int> &ends_1, list<int> &starts_2, list<int> &ends_2){
  //Combine the ancestries of the two nodes undergoing coalescence
	list<int> tempStarts;
	list<int> tempEnds;
	list<int>::iterator itStart1=starts_1.begin(), itStart2=starts_2.begin(), itEnd1=ends_1.begin(), itEnd2=ends_2.begin();
	while((itStart1 != starts_1.end()) || (itStart2 != starts_2.end())){
		if (itStart1 == starts_1.end()){//Finished searching first list
			add_interval_to_interval_list(tempStarts, tempEnds, *itStart2, *itEnd2);
			++itStart2;
			++itEnd2;
		}else if (itStart2 == starts_2.end()){//Finished searching second list
			add_interval_to_interval_list(tempStarts, tempEnds, *itStart1, *itEnd1);
			++itStart1;
			++itEnd1;
		}else{//Find which interval begins first and add to temporary list
			if (*itStart1 <= *itStart2){
				add_interval_to_interval_list(tempStarts, tempEnds, *itStart1, *itEnd1);
				++itStart1;
				++itEnd1;
			}else{
				add_interval_to_interval_list(tempStarts, tempEnds, *itStart2, *itEnd2);
				++itStart2;
				++itEnd2;
			}
		}
	}
  starts_1 = tempStarts;
  ends_1 = tempEnds;
}
#endif
