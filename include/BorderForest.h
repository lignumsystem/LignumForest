/// \file BorderForest.h
/// \brief Border forest radiation extinction 
#ifndef BORDERFOREST_H
#define BORDERFOREST_H
#include <iterator>
#include <Lignum.h>
using namespace Lignum;

namespace LignumForest{
  ///\brief Border forest radiation extinction
  ///
  ///BorderForest realizes the effect of forest outside the stand on radiation
  ////calculation. At the moment "green box" with Lambert Beer extinction
  ///in it is assumed.
  ///
  ///The corners of the stand are left (= minx, miny) and right (maxx, maxy), they
  ///correspond to lower left and upper right of voxelspace but are not necessarily
  ///the same since borders of the stand and borders of the voxelspace may not
  ///be the same.
  class BorderForest{
  public:
    BorderForest() : k_e(0.14){}
    BorderForest(const Point l, const Point r) : k_e(0.14), corner_l(l),
						 corner_r(r) {}
    void setH(LGMdouble h) {H = h;}
    void setHcb(LGMdouble hcb) {Hcb = hcb;}
    void setLAI(LGMdouble l) {LAI = l;}
    void setKExt(LGMdouble ke) {k_e = ke;}
    void setCornerL(const Point l) {corner_l = l;}
    void setCornerR(const Point r) {corner_r = r;}
    LGMdouble getH() {return H;}
    LGMdouble getHcb() {return Hcb;}
    LGMdouble getLAI() {return LAI;}
    LGMdouble getKExt() {return k_e;}
    Point& getCornerL() {return corner_l;}
    Point& getConrerR() {return corner_r;}
    ///\brief Extinction caused by the homogenous conifer border forest
    ///
    ///\pre  |\p dir| = 1
    ///\pre \p start Start point inside the voxel space
    ///\param start Start point of the light beam
    ///\param dir Direction of the light beam
    ///\param k_conifer The extinction coeffient value for homogenousus conifer border forest 
    ///\retval ext The extinction caused by the border stand.
    ///
    ///\par Border forest extinction 
    /// -# Calculate the  points where  the light beam  exits the  voxel space
    ///  - Check all six voxel space slabs for intersection to the direction \p dir and return the closest positive distance.
    ///   - Select the shortest distance \p tdist to the voxel space boundary.
    ///  - Calculate the \p exit point through one voxel space slab from the \p start point.
    ///   - The ray is represented as \f$ \mathrm{d_0}+\mathrm{tdist} \times \mathrm{dir}\f$.
    /// -# Use voxelspace::NearbyShading() to calculate and return the extinction value.
    ///  - \p BorderForest::H and \p BorderForest::Hcb define the border forest layer height.
    ///  - \p BorderForest::LAI and \p BorderForest::k_conifer specify border forest optical properties.
    ///  - Border forest is coniferous, LAI for deciduous trees and \p BorderForest::k_deciduous are 0.0.
    ///  .
    ///\snippet{lineno} borderforest.cc DEXT
    ///\sa voxelspace::NearbyShading()
    ///\par The C++ standard library: bind2nd deprecated 
    ///Replaced  \p std::bind2nd with \p std::bind.
    ///\code{.cpp}
    /// //Old deprecated version with bind2nd
    /// vector<double>::iterator it = find_if(v.begin(),v.end(),
    ///                                       bind2nd(greater_equal<double>(),0.0));
    /// //New version with bind
    /// vector<double>::iterator it = find_if(v.begin(),v.end(),
    ///                                       bind(greater_equal<double>(),std::placeholders::_1,0.0));
    ///\endcode
    ///The use of \p std::bind may not be obvious. \p std::bind creates function object with some number of parameters.
    ///The \p std::placeholders  denote how the actual parameters map to formal parameters in a function call.
    ///For example, the three function calls are equivalent:
    ///\code{.cpp}
    ///  fbind = bind(greater_equal<double>(),std::placeholders::_1,0.0)
    ///  fbind2nd = bind2nd(greater_equal<double>(),0.0));
    ///  fbind(5.0)
    ///  fbind2nd(5.0)
    ///  greater_equal<double>()(5.0,0.0)
    ///\endcode
    ///\note Modern C++ guidelines recommend C++ lambda expressions over std::bind for readability. 
    LGMdouble getBorderForestExtinction(const Point& start, 
					const PositionVector& dir, LGMdouble k_conifer);
  private:
    LGMdouble k_e; ///< Conifer extinction coefficent (\f$k\f$} value 
    LGMdouble H; ///< Crown height 
    LGMdouble Hcb; ///< Crown base height
    LGMdouble LAI; ///< Leaf area index
    Point corner_l, corner_r;  ///<corners of the stand that is bordered
  };
}//End namespace LignumForest
#endif
