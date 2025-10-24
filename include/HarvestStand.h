/// \file HarvestStand.h
/// \brief Self thinning in forest stand
#ifndef HARVEST_STAND_H
#define HARVEST_STAND_H

#include <mathsym.h>
#include <vector>
#include <Lignum.h>
#include <ScotsPine.h>
 
using namespace std;
using namespace cxxadt;

namespace LignumForest{
  ///\bref Remove a location with probability \p r to
  ///
  ///Remove a location with probability 'r' to get the wanted stand density
  ///\param v Boolean vector to denote removal
  ///\param r Removal probabilty
  ///\return Number of removals
  ///\pre \p v initialized with false
  ///\post \p v[i] = true if removal \p v[i] = false otherwise 
  int HarvestStand(vector<bool>& v, double r);
  ///\brief Remove trees in a gap
  ///
  ///Remove all trees that  are within the given  area of the tree located at point \p p.
  ///\param v Vector of tree locations
  ///\param rpos  Boolean vector denoting if a position is to be removed
  ///\param p position of the target tree
  ///\param gap_radius radius of the gap round the target tree
  ///\param verbose Verbose output
  ///\return Number of removed trees
  int ClearGap(vector<pair<double,double> >& v,  vector<bool>& rpos, const Point& p, double gap_radius, bool verbose);

  double SelfThinning(double Dbase, double Dbh, double N);

  double fN12(double d);

  ///\brief Remove trees.
  ///
  ///Remove and delete locations, trees and lsystems from the positions denoted 
  ///by the the vector \p rpos.
  ///\param rpos Boolean vector, *true* remove tree *false* maintain the tree
  ///\param locations Vector of tree locations
  ///\param vtree Vector of trees
  ///\param vlsystem Vector of L-systems for trees.
  ///\return \f$|locations|\f$
  ///\pre  \f$|rpos| = |locations| = |vtree| = |vlsystems|\f$
  ///\post \f$|rpos| = |locations| = |vtree| = |vlsystems|\f$
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
