#ifndef TERMINATE_ESCAPED_BUDS_H
#define TERMINATE_ESCAPED_BUDS_H
#include <RectangularCuboid.h>
#include <Lignum.h>
/// \file TerminateEscapedBuds.h
/// \brief Terminate buds not inside voxel space
namespace LignumForest{
  ///\brief Terminate tree buds that are not inside VoxelSpace.
  ///
  ///This functor can be used with ForEach algorithm to delete buds grown outside the VoxelSpace.
  ///\tparam TS TreeSegment
  ///\tparam BUD Bud
  template<class TS,class BUD> 
  class TerminateEscapedBuds{
  public:
    ///\brief Delegate constructor
    ///
    ///Create empty TerminateEscapedBuds::vs_orig and TerminateEscapedBuds::cs_orig
    TerminateEscapedBuds():TerminateEscapedBuds(Point(0,0,0),Point(0,0,0),0,0){}
    ///\brief Constructor with VoxelSpace corner points and border stand widths.
    ///
    ///Set TerminateEscapedBuds::vs_orig according to \p vs_ll and \p vs_ur.
    ///Set TerminateEscapedBuds::cs_orig corner points by adding \p x_border and \p y_border to \p vs_ll
    ///and subtracting them from \p vs_ur.
    ///\param vs_ll VoxelSpace lower left point, origo
    ///\param vs_ur VoxelSpace upper right point, diagionally opposite to \p vs_lower_right
    ///\param x_border Width of the border forest along x-axis, exterior portion of the center stand
    ///\param y_border Width of the border forest along y-axis, exterior portion of the center stand
    TerminateEscapedBuds(const Point& vs_ll, const Point& vs_ur, double x_border, double y_border);
    ///\brief Resize with VoxelSpace corner points and border stand widths.
    ///
    ///Set TerminateEscapedBuds::vs_orig  according to \p vs_ll and \p vs_ur.
    ///Set TerminateEscapedBuds::cs_orig corner points by adding \p x_border and \p y_border to \p vs_ll
    ///and subtracting them from \p vs_ur.
    ///\param vs_ll VoxelSpace lower left point, origo
    ///\param vs_ur VoxelSpace upper right point, diagionally opposite to \p vs_lower_right
    ///\param x_border Width of the border forest along x-axis, exterior portion of the center stand
    ///\param y_border Width of the border forest along y-axis, exterior portion of the center stand
    ///\retval self Resized TerminateEscapedBuds
    ///\note The same algorithm as in the constructor TerminateEscapedBuds(const Point&, const Point&, double, double);
    ///\sa TerminateEscapedBuds::TerminateEscapedBuds(const Point&, const Point&, double, double)
    TerminateEscapedBuds& resize(const Point& vs_ll, const Point& vs_ur, double x_border, double y_border);
    ///\brief Terminate a bud if it is outside VoxelSpace TerminateEscapedBuds::vs_orig
    ///\param tc A tree compartment
    ///\pre The tree compartment \p tc is a Bud
    ///\post The Bud state is DEAD if it is outside VoxelSpace
    ///\retval tc The tree compartment
    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;
    ///\brief Check if a point is inside center stand TerminateEscapedBuds::cs_orig 
    ///\param p The point to examine
    ///\return True if \p p inside the center stand, False othewise
    bool insideCenterStand(const Point& p)const{return cs_orig.insideVolume(p);}
    ///\brief Check if a point is inside VoxelSpace TerminateEscapedBuds::vs_orig
    ///\param p The point to examine
    ///\return True if \p p inside the VoxelSpace, False othewise
    bool insideVoxelSpace(const Point& p)const{return vs_orig.insideVolume(p);}
  private:
    RectangularCuboid vs_orig;///< Original VoxelSpace perimeter
    RectangularCuboid cs_orig;///< Original center stand perimeter
  };

  template<class TS,class BUD>
  TerminateEscapedBuds<TS,BUD>::TerminateEscapedBuds(const Point& vs_ll, const Point& vs_ur, double x_border, double y_border)
  {
    vs_orig.resize(vs_ll,vs_ur);
    double x_ll = vs_ll.getX();
    double y_ll = vs_ll.getY();
    //Add the exterior portion to get lower left point for center stand
    Point cs_ll = Point(x_ll+x_border,y_ll+y_border,vs_ll.getZ());
    double x_ur = vs_ur.getX();
    double y_ur = vs_ur.getY();
    //Subtract the exterior portion to get upper right point for center stand
    Point cs_ur = Point(x_ur-x_border,y_ur-y_border,vs_ur.getZ());
    cs_orig.resize(cs_ll,cs_ur);
  }

  template<class TS,class BUD>
  TerminateEscapedBuds<TS,BUD>& TerminateEscapedBuds<TS,BUD>::resize(const Point& vs_ll, const Point& vs_ur, double x_border, double y_border)
  {
    vs_orig.resize(vs_ll,vs_ur);
    double x_ll = vs_ll.getX();
    double y_ll = vs_ll.getY();
    //Add the exterior portion to get lower left point for center stand
    Point cs_ll = Point(x_ll+x_border,y_ll+y_border,vs_ll.getZ());
    double x_ur = vs_ur.getX();
    double y_ur = vs_ur.getY();
    //Subtract the exterior portion to get upper right point for center stand
    Point cs_ur = Point(x_ur-x_border,y_ur-y_border,vs_ur.getZ());
    cs_orig.resize(cs_ll,cs_ur);
    return *this;
  }
    
  template<class TS,class BUD>
  TreeCompartment<TS,BUD>* TerminateEscapedBuds<TS,BUD>::operator()(TreeCompartment<TS,BUD>*tc)const
  {
    if (BUD* b = dynamic_cast<BUD*>(tc)){
      Point p = GetPoint(*b);
      if (!insideVoxelSpace(p)){
	SetValue(*b,LGAstate,DEAD);
      }
    }
    return tc;
  }
}// namespace LignumForest
#endif
