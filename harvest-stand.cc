#include <HarvestStand.h>
namespace LignumForest{
//Stand harvest: remove a location with probability 'r' to
//get the wanted stand density 
//rpos boolean vector denoting if a position is to be removed
//r probability to be removed
//return: 
//rpos vector of boolean values denoting if the position is to
//be removed or not
int  HarvestStand(vector<bool>& rpos, double r)
{
  int r_seed = 1;
  int removed = 0;
  for (unsigned int i = 0; i < rpos.size(); i++){
    double p = ran3(&r_seed);
    if (p <= r){
      rpos[i] = true;
      removed++;
    }
  }
  //the final density
  return removed;
}

//Clear gap:  remove all trees that  are within the given  area of the
//tree located at point p.
//v vector of tree locations
//rpos  boolean vector denoting if a position is to be removed
//p position of the target tree
//gap_radius radius of the gap round the target tree
int ClearGap(vector<pair<double,double> >& v, vector<bool>& rpos, const Point& p, 
	     double gap_radius, bool verbose)
{
  int removed = 0;
  vector<pair<double,double> >::iterator current = v.begin();
  for (unsigned int i = 0; i < rpos.size(); i++){
    pair<double,double>& pair = v[i];
    const Point l = Point(pair.first,pair.second,0.0);
    if ((l || p) <= gap_radius){
      if (verbose){
	cout << "ClearGap, removing tree: " << "  Radius: " <<  gap_radius << " Distance: " 
	     << (l||p) << " Point: " << l <<endl;
      }
      rpos[i] = true;
      removed++;
    }
  }
  return removed;
} 
}

