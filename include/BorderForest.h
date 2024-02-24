#ifndef BORDERFOREST_H
#define BORDERFOREST_H
#include <iterator>
#include <Lignum.h>
using namespace Lignum;

namespace LignumForest{
  ///\brief Border forest extinction
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
    ///\brief Extinction of homogenous conifer border forest
    ///\param start   Start point of the light beam
    ///\param dir  Direction of the light beam, \note |dir| == 1
    ///\param k_conifer The K value for homogenousus conifer border forest 
    ///\retval tau The extinction caused by the border stand
    ///Calculate the  point where  the light beam  exits the  voxel space
    ///(there  must be  one). NearByShading  then returns  the extinction
    ///coeffcient. \sa NearByShading
    LGMdouble getBorderForestExtinction(const Point& start, 
					const PositionVector& dir, LGMdouble k_conifer);
  private:
    LGMdouble k_e; ///< Conifer K value 
    LGMdouble H; ///< Crown height 
    LGMdouble Hcb; ///< Crown base height
    LGMdouble LAI; ///< Leaf area index
    Point corner_l, corner_r;  ///<corners of the stand that is bordered
  };
}//End namespace LignumForest
#endif
