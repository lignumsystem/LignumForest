#include <TreeLocations.h>
///\file generate-tree-locations.cc
///\brief Tree locations  Nonstationary Poisson process.
extern int ran3_seed;   //is a global variable
namespace LignumForest{
/// Generate tree locations: Nonstationary Poisson process: Stoyan, Kendall & Mecke, p. 52 -54
/// \param nTrees    number of trees (input, output)
/// \param corner1X  lower left  corner  of forest,  X, Y  (input)
/// \param corner1Y  lower left  corner  of forest,  X, Y  (input)
/// \param corner2X  upper  right  corner  (input)  
/// \param corner2Y  upper  right  corner  (input)  
/// \param rmin      minimum distance between  the trees  (input) 
/// \param gap       representation of a circular gap: ForestGap is simply pair<pair<x,y>,r> 
///          where pair<x,y> is the position and r the radius of the circular gap
/// \param v         vector  containing accepted x,y tree location  coordinates (output)  
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

  ///First: generate all coordinates 
  for(int i = 0; i < nTrees; i++)   {
    xVal[i] = corner1X + xDist * ran3( &ran3_seed );
    yVal[i] = corner1Y + yDist * ran3( &ran3_seed );
  }


  ///Second (hard core): delete one tree from each pair of trees nearer than
  ///a minimum distance. 
  vector<bool> iDeleted(nTrees);
  for(int i = 0; i < nTrees; i++)
    iDeleted[i] = false;
  
  ///Third: pairwise comparison (almost), if a  problem (nTrees > 1.0e+04), 
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
  
  ///Fourth: after deletion copy the x- and y- coordinates of not deleted
  ///trees to items 0 ,.., new-value-of-nTrees of xVal and yVal
  for(int i = 0; i < nTrees; i++)
    if(!iDeleted[i])        {
      v.insert(v.end(),pair<double,double>(xVal[i],yVal[i]));
    }
  ///Finally: return number of positions actually created (the first argument)
  nTrees = v.size();
}

}
