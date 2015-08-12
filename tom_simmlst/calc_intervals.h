void calc_intervals(const vector<bool> &AncMat, vector<int> &intervalStarts, vector<int> &intervalEnds, const vector<int> &blockStarts, const vector<int> &blockEnds, const int b){
  //Calculate the intervals of ancestral material for the given node
  if (b == 1){
    //Only one block of ancestral material
    if (AncMat[blockStarts[0]] == true){
      intervalStarts.push_back(blockStarts[0]);
      //Check for interval of length 1
      if (AncMat[blockStarts[0]+1] == false) intervalEnds.push_back(blockStarts[0]);
    }
    for (int a=blockStarts[0]+1;a<blockEnds[0]-1;++a){
      if ((AncMat[a] == true) && (AncMat[a-1] == false)){
        intervalStarts.push_back(a);
        //Check for interval of length 1
        if (AncMat[a+1] == false) intervalEnds.push_back(a);
      }else if ((AncMat[a] == true) && (AncMat[a+1] == false)) intervalEnds.push_back(a);
    }
    //Check last element
    if (AncMat[blockEnds[0]-1] == true){
      intervalEnds.push_back(blockEnds[0]-1);
      //Check for interval of length 1
      if (AncMat[blockEnds[0]-2] == false) intervalStarts.push_back(blockEnds[0]-1);
    }

  }else{
    //More than one ancestral block
    //First ancestral block
    if (AncMat[blockStarts[0]] == true){
      intervalStarts.push_back(blockStarts[blockStarts[0]]);
      //Check for interval of length 1
      if (AncMat[blockStarts[0]+1] == false) intervalEnds.push_back(blockStarts[0]);
    }
    for (int a=blockStarts[0]+1;a<blockEnds[0];++a){
      if ((AncMat[a] == true) && (AncMat[a-1] == false)){
        intervalStarts.push_back(a);
        //Check for interval of length 1
        if (AncMat[a+1] == false) intervalEnds.push_back(a);
      }else if ((AncMat[a] == true) && (AncMat[a+1] == false)) intervalEnds.push_back(a);
    }

    //Middle ancestral blocks
    for (int m=1;m<b-1;++m){
      //Find start and end points
      for (int a=blockStarts[m];a<blockEnds[m];++a){
        if ((AncMat[a] == true) && (AncMat[a-1] == false)){
          intervalStarts.push_back(a);
          //Check for interval of length 1
          if (AncMat[a+1] == false) intervalEnds.push_back(a);
        }else if ((AncMat[a] == true) && (AncMat[a+1] == false)) intervalEnds.push_back(a);
      }
    }

    //Last ancestral block
    for (int a=blockStarts[b-1];a<blockEnds[b-1]-1;++a){
      if ((AncMat[a] == true) && (AncMat[a-1] == false)){
        intervalStarts.push_back(a);
        //Check for interval of length 1
        if (AncMat[a+1] == false) intervalEnds.push_back(a);
      }else if ((AncMat[a] == true) && (AncMat[a+1] == false)) intervalEnds.push_back(a);
    }
    //Check last element of block
    if (AncMat[blockEnds[b-1]-1] == true){
      intervalEnds.push_back(blockEnds[b-1]-1);
      //Check for interval of length 1
      if (AncMat[blockEnds[b-1]-2] == false) intervalStarts.push_back(blockEnds[b-1]-1);
    }
  }
}
