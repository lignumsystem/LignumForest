#ifndef HARVEST_STAND_H
#define HARVEST_STAND_H

#include <mathsym.h>
#include <vector>
#include <Lignum.h>
#include <ScotsPine.h>
 
using namespace std;
using namespace cxxadt;

namespace LignumForest{ 
  int HarvestStand(vector<bool>& v, double r);
  int ClearGap(vector<pair<double,double> >& v,  vector<bool>& rpos, const Point& p, double gap_radius, bool verbose);

  double SelfThinning(double Dbase, double Dbh, double N);

  double fN12(double d);

  ///Remove and delete locations, trees and lsystems from the positions denoted 
  ///by the the vector rpos
  template <class TREE, class LSYSTEM>
  unsigned int RemoveTrees(vector<bool>& rpos,vector<pair<double,double> >& locations, 
			   vector<TREE*>& vtree, vector<LSYSTEM*>& vlsystem)
  {
    typename vector<bool>::iterator rpos_it = rpos.begin();
    typename vector<pair<double,double> >::iterator locations_it = locations.begin(); 
    typename vector<TREE*>::iterator vtree_it  = vtree.begin();
    typename vector<LSYSTEM*>::iterator vlsystem_it = vlsystem.begin();
    //A single tree, a copy of it is inserted into voxel space
    //Remove positions from the locations vector 
    if (vtree.size() == 1){
      while (rpos_it != rpos.end()){
	if (*rpos_it == true){
	  rpos.erase(rpos_it);
	  locations.erase(locations_it);
	}
	else{
	  advance(rpos_it,1);
	  advance(locations_it,1);
	}
      }
    }
    //Many individual trees, remove locations, trees and L-systems
    //Delete trees and L-systems
    else{
      while(rpos_it != rpos.end()){
	if (*rpos_it == true){
	  rpos.erase(rpos_it);
	  locations.erase(locations_it);
	  delete *vtree_it;
	  vtree.erase(vtree_it);
	  delete *vlsystem_it;
	  vlsystem.erase(vlsystem_it);
	}
	else{
	  advance(rpos_it,1);
	  advance(locations_it,1);
	  advance(vtree_it,1);
	  advance(vlsystem_it,1);
	}
      }
    }
    return locations.size();
  }

}//End namespace LignumForest
      
#endif
