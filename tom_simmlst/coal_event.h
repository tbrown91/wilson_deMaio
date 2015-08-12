#ifndef COAL_EVENT
#define COAL_EVENT
// add an interval to a list, or join it to the last interval of the list if they overlap
void add_interval_to_interval_list(vector<int> &starts, vector<int> &ends, int start, int end){
	int index=starts.size()-1;
  if (index == -1){
    starts.push_back(start);
    ends.push_back(end);
  }else{
  	if (start<=(ends[index]+1)){
  		if (end>ends[index]){
  			ends[index]=end;
  		}
  	}else{
  		starts.push_back(start);
  		ends.push_back(end);
  	}
  }
}
void combine_ancestries(vector<int> &starts_1, vector<int> &ends_1, vector<int> &starts_2, vector<int> &ends_2){
  //Combine the ancestries of the two nodes undergoing coalescence
  vector<int> tempStarts;
  vector<int> tempEnds;
  int index_1=0, index_2=0;
  int b_1 = starts_1.size();
  int b_2 = starts_2.size();
  while (true){
    if ((index_1 == b_1) && (index_2 == b_2)) break;
    else if (index_1 == b_1){
      add_interval_to_interval_list(tempStarts, tempEnds, starts_2[index_2], ends_2[index_2]);
      ++index_2;
    }else if (index_2 == b_2){
      add_interval_to_interval_list(tempStarts, tempEnds, starts_1[index_1], ends_1[index_1]);
      ++index_1;
    }else{
      //Find which of the current intervals starts first
      if (starts_1[index_1] <= starts_2[index_2]){
      	add_interval_to_interval_list(tempStarts, tempEnds, starts_1[index_1], ends_1[index_1]);
      	++index_1;
      }else{
      	add_interval_to_interval_list(tempStarts, tempEnds, starts_2[index_2], ends_2[index_2]);
      	++index_2;
      }
    }
  }
  starts_1 = tempStarts;
  ends_1 = tempEnds;
}
#endif
