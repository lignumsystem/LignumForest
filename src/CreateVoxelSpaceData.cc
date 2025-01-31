#include <CreateVoxelSpaceData.h>
namespace LignumForest{
  void CreateVoxelSpaceData::insertData(const VoxelSpace& vs,int year,int interval){
    if (year % interval == 0){
      Point lower_left = vs.getLowerLeftCorner();
      Point upper_right = vs.getUpperRightCorner();
      double llx,lly,llz,urx,ury,urz=0.0;
      vector<double> v(7,0);
      llx = lower_left.getX();
      lly = lower_left.getY();
      llz = lower_left.getZ();
      urx = upper_right.getX();
      ury = upper_right.getY();
      urz = upper_right.getZ();
      v[0] = static_cast<double>(year);
      v[1] = llx;
      v[2] = lly;
      v[3] = llz;
      v[4] = urx;
      v[5] = ury;
      v[6] = urz;
      vspacedata.append(v);
      current_row = current_row + 1;
    }
  }
}
  
