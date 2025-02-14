#include <TreeLocations.h>
#include <LignumForestGlobals.h>
using namespace LignumForest;
///\file generate-tree-locations.cc
//extern int LignumForest::ran3_seed;   //is a global variable
namespace LignumForest{
  ///\brief Tree locations  Nonstationary Poisson process.
  void GenerateLocations(int& nTrees, double corner1X, double corner1Y, double corner2X, 
			 double corner2Y, double rmin, 
			 const ForestGap& gap,
			 vector<pair<double,double> >& v)
  { 


    double xDist = corner2X - corner1X;
    double yDist = corner2Y - corner1Y;
  
    //xVal     array containing x coordinates of tree locations 
    //yVal     array containing y coordinates of tree locations 
    vector<double> xVal(nTrees);
    vector<double> yVal(nTrees);

    ///\par First: generate all coordinates 
    for(int i = 0; i < nTrees; i++)   {
      xVal[i] = corner1X + xDist * ran3( &LignumForest::ran3_seed );
      yVal[i] = corner1Y + yDist * ran3( &LignumForest::ran3_seed );
    }


    ///\par Second (hard core): delete one tree from each pair of trees nearer than
    ///a minimum distance. 
    vector<bool> iDeleted(nTrees);
    for(int i = 0; i < nTrees; i++)
      iDeleted[i] = false;
  
    ///\par Third: pairwise comparison (almost), if a  problem (nTrees > 1.0e+04), 
    ///redesign the  generation of the  coordinates so that they  will be
    ///inserted into the vector only if the distance is acceptable.
    for(int i = 0; i < nTrees - 1; i++){
      if(!iDeleted[i])
	for(int j = i + 1; j < nTrees; j++){
	  if(!iDeleted[j]){
	    //check the distance
	    if(pow(xVal[i]-xVal[j], 2.0) + pow(yVal[i]-yVal[j], 2.0)
	       < rmin * rmin){
	      iDeleted[j] = true;
	    }
	    //check the gap
	    if (pow(gap.first.first-xVal[j], 2.0) + pow(gap.first.second-yVal[j],2.0)
		< gap.second*gap.second){
	      iDeleted[j] = true; 
	    }
	  }
	}//for (int j
    }//for (int i 
  
    ///\par Fourth: after deletion copy the x- and y- coordinates of not deleted
    ///trees to items 0 ,.., new-value-of-nTrees of xVal and yVal
    for(int i = 0; i < nTrees; i++)
      if(!iDeleted[i])        {
	v.insert(v.end(),pair<double,double>(xVal[i],yVal[i]));
      }
    ///\par Finally: return number of positions actually created (the first argument)
    nTrees = v.size();
  }
}//End namespace LignumForest
