#ifndef UPDATE_MRCA
#define UPDATE_MRCA
//#include "modifyMRCA.h"
//Update MRCA structure providing intervals that are shared by the two coalescing lineages
void update_MRCA(ARG::MRCA &M, vector<int> &starts1, vector<int> &ends1, vector<int> &starts2, vector<int> &ends2){
  int index1=0, index2=0;
  int b1 = starts1.size(), b2 = starts2.size();
  (M).itStart=(M).starts.begin();
  (M).itEnd=(M).ends.begin();
  (M).itValue=(M).values.begin();
  int currentStart1=0, currentStart2=0;
  if ((index1 != b1) && (index2 != b2)) currentStart1=starts1[index1], currentStart2=starts2[index2];
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
    	//modifyMRCA(M, currentStart1, ends1[index1]);
    	currentStart2=ends1[index1]+1;
    	++index1;
    	if (index1<b1) currentStart1=starts1[index1];
    }else if ((currentStart1 == currentStart2) && (ends2[index2]<ends1[index1])){ //overlap, same start
    	//modifyMRCA(M, currentStart2, ends2[index2]);
    	currentStart1=ends2[index2]+1;
    	++index2;
    	if (index2<b2) currentStart2=starts2[index2];
    }else if ((currentStart1 == currentStart2) && (ends2[index2]==ends1[index1])){ //overlap, same start, same end
    	//modifyMRCA(M, currentStart2, ends2[index2]);
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
#endif
