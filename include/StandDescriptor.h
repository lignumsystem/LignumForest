#ifndef STAND_DESCRIPTOR_H
#define STAND_DESCRIPTOR_H
#include <LGMSymbols.h>
using namespace Lignum;
/// \file StandDescriptor.h

/// \brief StandDescriptor collects and evaluates stand level quantities.
///
/// It evalueates the the stand level values from a vector of trees (input).
/// Query of values is now accomplished with methods (getXXX()). This may later
/// be changed to GetValue construct.
/// \tparam TREE Lignum tree
template <class TREE>
class StandDescriptor{
public:
  ///This constructor sets hard coded area and corner points  
  StandDescriptor() :
    dbh_mean(0.0), hdom(0.0), stemVolume(0.0),
    standBasalArea(0.0), meanHeight(0.0), basalAreaAtCrownBase(0.0), minDBH(R_HUGE),
    maxDBH(-R_HUGE), noTrees(0.0), age(0), minH(R_HUGE), maxH(-R_HUGE), LAI(0.0)
  { corner_ll=Point(0.0,0.0,0.0);
    corner_ur= Point(10.0,10.0,0.0);
    area = 10.0*10.0/10000.0;
  }

  ///This constructor uses given corner points for area (in voxel space)
  ///\sa VoxelSpace 
  StandDescriptor(const Point& ll, const Point& ur):
    corner_ll(ll), corner_ur(ll),
    dbh_mean(0.0), hdom(0.0), stemVolume(0.0),
    standBasalArea(0.0), meanHeight(0.0), basalAreaAtCrownBase(0.0), minDBH(R_HUGE),
    maxDBH(-R_HUGE), minH(R_HUGE), maxH(-R_HUGE), age(0), LAI(0.0), noTrees(0.0)
  {area = (ur.getX()-ll.getX())*(ur.getY()-ll.getY());} //assume rectangle

  LGMdouble getMeanDbh()  { return dbh_mean; }
  LGMdouble getMeanDbase() {return dbase_mean;}
  LGMdouble getHDom()   { return hdom; }
  LGMdouble getMeanHeight() { return meanHeight; }
  LGMdouble getStemVolume()	{ return stemVolume; }
  LGMdouble getStandBasalArea() { return standBasalArea; }
  LGMdouble getMeanDensity() {
    if(area > 0.0)
      return (double)noTrees/area;
    else
      return 0.0;
  }
  int getNoTrees() {return noTrees;}
  LGMdouble getMinCrownLimit(){ return minCrownLimit; }
  LGMdouble getMaxCrownLimit(){ return maxCrownLimit; }
  LGMdouble getMeanCrownLimit(){ return meanCrownLimit; }
  LGMdouble getArea() { return area; }
  int getAge()	{ return age; }
  LGMdouble getBasalAreaAtCrownBase() { return basalAreaAtCrownBase; }
  Point getCornerLl( ) {return corner_ll;}
  Point getCornerUr() {return corner_ur;};
  LGMdouble getMinDBH() { return minDBH; }
  LGMdouble getMaxDBH() { return maxDBH; }
  LGMdouble getMinH() {return minH;}
  LGMdouble getMaxH() {return maxH;}
  LGMdouble getMinDbase() {return minDbase; }
  LGMdouble getMaxDbase() {return maxDbase; }
  LGMdouble getLAI() {return LAI;}
  LGMdouble getWfMass() {return Wf;}
  /// \brief Evaluate stand level variables
  /// \param vt Vector of trees (from StandDescriptor)
  /// \param loc Vector of tree locations (from StandDescriptor)
  /// \tparam TREE Lignum tree
  /// \sa GrowthLoop::evaluateStandVariables
  /// \sa GrowthLoop::vtree GrowthLoop::locations
  void evaluateStandVariables(vector<TREE*> vt, vector<pair<double,double> >& loc);
  void setAge(int a) {age = a;}
  void setLlCorner(const Point& c) {corner_ll = c;}
  void setUrCorner(const Point& c) {corner_ur = c;}
  void evaluateArea() {area = (corner_ur.getX()-corner_ll.getX())*(corner_ur.getY()-corner_ll.getY());}
  /// \brief Output stand data to a file
  /// \param stand_output Output file
  /// \pre Stand variables have been evaluated
  /// \sa evaluateStandVariables
  void writeOutput(ofstream& stand_output);
  bool inPlot(const Point& p) {
    if(p.getX() >= corner_ll.getX() &&
       p.getY() >= corner_ll.getY() &&
       p.getX() <= corner_ur.getX() &&
       p.getY() <= corner_ur.getY())
      return true;
    else
      return false;
  }

   

private:
  Point corner_ll, corner_ur; ///< Lower left and upper right corners
  LGMdouble dbh_mean, hdom, stemVolume;
  LGMdouble minCrownLimit, standBasalArea;
  LGMdouble maxCrownLimit, meanCrownLimit;
  LGMdouble meanHeight;
  LGMdouble basalAreaAtCrownBase;
  LGMdouble minDBH, maxDBH;
  int  noTrees, age;
  LGMdouble minH, maxH, LAI;
  LGMdouble area;  ///< Area in hectars
  LGMdouble dbase_mean;
  LGMdouble minDbase, maxDbase;
  LGMdouble Wf;
};

#include <StandDescriptorI.h>

    

#endif
