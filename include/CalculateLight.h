/// \file CalculateLight.h
/// \brief Radiation calculations
///
/// Two functors for Scots pine
/// \arg \c EvaluateRadiationForCfTreeSegmentInVoxelSpace
/// \arg \c AccumulateOpticalDepth
#ifndef CALCULATE_LIGHT_H
#define CALCULATE_LIGHT_H
#include <Lignum.h>
#include <ScotsPine.h>
#include <CopyDumpCfTree.h>
#include <VoxelSpace.h>
#include <SomeFunctors.h>
#include <vector>
#include <utility>
#include <BorderForest.h>
#include<limits>
using namespace std;

using namespace Lignum;
using namespace sky;


#define HIT_THE_FOLIAGE 1
#define NO_HIT 0
#define HIT_THE_WOOD -1

//===============================================================================
// VoxelSpace calculation
namespace LignumForest{
  
  template <class TS, class BUD>
  class EvaluateRadiationForCfTreeSegmentInVoxelSpace {
  public:
    EvaluateRadiationForCfTreeSegmentInVoxelSpace(const ParametricCurve& k) : K(k),
									      evaluate_border(false) {}
    EvaluateRadiationForCfTreeSegmentInVoxelSpace(const ParametricCurve& k,
						  VoxelSpace* vs,BorderForest* bf, bool border,
						  LGMdouble kbc, bool wd, bool ps):
      K(k),voxel_space(vs), border_forest(bf), evaluate_border(border),
      k_border_conifer(kbc), wood(wd), pairwise_self(ps) {}

    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;

  private:
    const ParametricCurve& K;
    VoxelSpace* voxel_space;
    BorderForest* border_forest;
    bool evaluate_border;
    LGMdouble k_border_conifer;
    bool wood;         //If the needleless woody parts are included into calculation
    bool pairwise_self; //If self-shading is evaluated by ray casting
  };

  class AccumulateOpticalDepth{
  public:
    AccumulateOpticalDepth(LGMdouble side, bool wd) :
      box_side_length(side), wood(wd)
    {
      box_volume = pow(box_side_length,3.0);
    }

    double operator()(double o_d,VoxelMovement& vm){
      if(vm.af > R_EPSILON || (wood && vm.wood_area > R_EPSILON)) {
	LGMdouble k = vm.STAR_mean;
	o_d += k * vm.af * vm.l / box_volume;
	if(wood) {
	  //Mean projection area of surface of a circular cylinder (excluding end disks)
	  // is 1/4 of its area
	  o_d += 0.25 * vm.wood_area * vm.l / box_volume;
	}
      }
      return o_d;
    }
  private:
    LGMdouble box_side_length;
    LGMdouble box_volume;
    bool wood;                       //  If woody parts are considered
  };

#undef HIT_THE_FOLIAGE
#undef NO_HIT
#undef HIT_THE_WOOD
}//End namespace LignumForest
#include <CalculateLightI.h>

#endif
