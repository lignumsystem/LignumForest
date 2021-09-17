#ifndef BORDERFOREST_H
#define BORDERFOREST_H
#include <iterator>
#include <Lignum.h>
using namespace Lignum;


//BorderForest realizes the effect of forest outside the stand on radiation
//calculation. At the moment "green box" with Lambert Beer extinction
//in it is assumed.

//The corners of the stand are left (= minx, miny) and right (maxx, maxy), they
//correspond to lower left and upper right of voxelspace but are not necessarily
//the same since borders of the stand and borders of the voxelspace may not
//be the same.
namespace Lignum{ 
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
    LGMdouble getBorderForestExtinction(const Point& start, 
					const PositionVector& dir, LGMdouble k_conifer);
 private:
    LGMdouble k_e;
    LGMdouble H;
    LGMdouble Hcb;
    LGMdouble LAI;
    Point corner_l, corner_r;   //corners of the stand that is bordered
};
}
#endif
