#include <TreeDataAfterGrowth.h>
#include <CreateVoxelSpaceData.h>
namespace LignumForest{
  void CreateVoxelSpaceData::insertData(const VoxelSpace& vs,int year,int interval){
    if (year % interval == 0){
      //Origo
      Point p0 = vs.getLowerLeftCorner();
      //Diaginally opposite point, clcokwise upper third point
      Point p6= vs.getUpperRightCorner();
      //Clockwise bottom third point
      Point p2(p6.getX(),p6.getY(),p0.getZ());
      //Clockwise bottom fourth point
      Point p3(p6.getX(),p0.getY(),p0.getZ());
      double llx,lly,llz,urx,ury,urz,p0p3,p2p3,p2p6,area=0.0;
      //The number of columns required  is the number of columns names
      vector<double> v(vs_sizes_column_names.size(),0);
      llx = p0.getX();
      lly = p0.getY();
      llz = p0.getZ();
      urx = p6.getX();
      ury = p6.getY();
      urz = p6.getZ();
      area = vs.getArea();
      //Width (x direction) of the voxel space
      p0p3 = p0||p3;
      //Length (y direction) of the voxels space
      p2p3 = p2||p3;
      //Height (z direction) of the voxel space
      p2p6 = p2||p6;
      v[0] = static_cast<double>(year);
      v[1] = llx;
      v[2] = lly;
      v[3] = llz;
      v[4] = urx;
      v[5] = ury;
      v[6] = urz;
      v[7] = area;
      v[8] = p0p3;
      v[9] = p2p3;
      v[10] = p2p6;
      v[11] = p0p3*p2p3;
      vspacedata.append(v);
      current_row = current_row + 1;
    }
  }
}
  
