#ifndef ARG_H
#define ARG_H
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <gsl/gsl_randist.h>
#include "data.h"
#include "popsize.h"
#include <string.h>
extern gsl_rng * rng;

using namespace std;

/**
    @brief This class represents a "small" ARG, that is one where all lines carry some ancestral material
*/
class Arg {

    public:
      Arg(int n,double rho,double rho_ext,double delta,double delta_ext,vector<int> blocks,vector<int> gaps,PopSize * popsize=NULL);///<Creates an ARG for n isolates, with recombination rate rho, tract length delta and block structure as given by the vector "blocks"
      Data * drawData(double theta,double theta_extMin, double theta_extMax);///<Create sequence data for the leaves of the ARG using theta/2 as mutation rate
      string extractCG();///<Extracts the clonal genealogy of the ARG
      string extractLT(int site);///<Extracts the local tree at the given site
      void outputDOT(ostream * out,bool am);///<Create a DOT description of the ARG and export it
      void outputLOCAL(ostream * out);///<Export the local trees
      void outputBREAKS(ostream * out);///<Write recombinant break intervals to log file
      static vector<int> makeBlocks(string str) {
        vector<int> v;
        int s=0;
        char * pch;
        pch = strtok ((char*)str.data(),",");
        v.push_back(s);
        while (pch!=NULL) {
            s+=atoi(pch);
            v.push_back(s);
            pch = strtok (NULL, ",");
          };
        return v;
      }
      static vector<int> makeGaps(string str) {
        vector<int> v;
        int s=0;
        char * pch;
        pch = strtok ((char*)str.data(),",");
        while (pch!=NULL) {
            s=atoi(pch);
            v.push_back(s);
            pch = strtok (NULL, ",");
          };
        return v;
      }
      // MRCA Keeps track along the genome of how many lineages (in values) have that interval (limited by starts and ends) in the ancestral material. Sites that reach MRCA will have a value of 1.
      struct MRCA {
          list<int> starts;
          list<int> ends;
          list<int> values;
          list<int>::iterator itStart;
          list<int>::iterator itEnd;
          list<int>::iterator itValue;
      };
    protected:
      void construct();///<Construction of the ARG, called by class constructor
      int n;///<Number of isolates
      double rho;///<Scaled recombination rate
      double rho_ext;//<External recombination rate
      double delta;///<Average tract length
      double delta_ext;///<Average tract length for external recombination
      vector<int> blocks;///<Structure of the observed data
      vector<int> gaps;///<Length of gap between each fragment of ancestral material
      vector<vector<int> > s;///<List of the nodes in the ARG, with s[.][0] and s[.][1] being the two children, s[.][2] and s[.][3] being the two parents, s[.][4] and s[.][5] being the start and end point of an import
      vector<double> ages;///<Ages of the nodes in the ARG
      vector<bool> clonal;///<Whether a node is part of the clonal genealogy or not
      vector<bool> changeLT;///<Indicates whether the local tree change at the sites

      // Iterators point at the particular interval that is considered at preset (useful for functions that explore it).
      string buildTree(int r);
      PopSize * popsize;///<Model for population size
  };
#endif
