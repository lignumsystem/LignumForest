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


// Functions in this file are to be used instead of radiation calculations in stl-lignum
// (stl-lignum/include/Shading.h * and ...I.h). Replace them eventually with this.
// _1 addition to the names are to distinquish the from functions of stl-lignum radiation
// calculations.

// EvaluateRadiationForCfTreeSegment_1 evaluates shading
// caused by all other segments on this conifer segment.  The shading caused by segments in the crown
// of tree itself is evaluated by ShadingEffectOfCfTreeSegment_1<TS,BUD> (it is the same as
// ShadingEffectOfCfTreeSegment in stl-lignum; the functor is duplicated here only for convenience).
// After that other trees are accounted for with voxelspace (voxel_space->getRoute() etc) and surrounding
// stand with border_forest->getBorderForestExtinction().
// This function does not evaluate the shading by surrounding trees (Qin_stand) and border forest separately 
// (slighthly faster);

//If the attributes voxel_space and border_forest are set (then
//voxel_space != NULL), the attenuation of the beam in the voxel_space
//and border_forest is taken into consideration.

template <class TS, class BUD>
class EvaluateRadiationForCfTreeSegment_1 {
public:
    EvaluateRadiationForCfTreeSegment_1(const ParametricCurve& k) : K(k),
        evaluate_voxel_and_border(false) {}
    EvaluateRadiationForCfTreeSegment_1(const ParametricCurve& k,
                                        VoxelSpace* vs,BorderForest* bf, LGMdouble gk):
        K(k),voxel_space(vs), border_forest(bf), evaluate_voxel_and_border(true),
        green_voxel_ext_coeff(gk) {}

    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;
private:
    const ParametricCurve& K;
    VoxelSpace* voxel_space;
    BorderForest* border_forest;
    bool evaluate_voxel_and_border;
    LGMdouble green_voxel_ext_coeff;
};


//This functor ShadingEffectOfCfTreeSegment<TS,BUD> evaluates shading caused
//by a conifer segment on this conifer segment (shaded_s)

template <class TS,class BUD>
class ShadingEffectOfCfTreeSegment_1 {
public:
    ShadingEffectOfCfTreeSegment_1(CfTreeSegment<TS,BUD>* ts, const ParametricCurve& K_in,
                                   vector<double>& sectors)
        :shaded_s(ts), K(K_in),S(sectors){}
    //ForEach functor to compute shadiness
    TreeCompartment<TS,BUD>*  operator()(TreeCompartment<TS,BUD>* tc)const;
    //Get vector for S (shadiness)
    vector<double>& getS(){return S;}
private:
    CfTreeSegment<TS,BUD>* shaded_s;
    //Avoid unnecessary constructor calls in generic algorithms
    const ParametricCurve& K;
    vector<double>& S;
};


//=======================================================================================================
//This version of radiation evaluates radiation conditions for subject tree by pairwise
// comparison

// EvaluateRadiationForCfTreeSegment_2 evaluates shading by all other segments by paiwise comparison (segments
// in own crown & other trees). It uses ShadingEffectOfCfTreeSegment_1<TS,BUD> to evaluate shading.
// This functor evaluates shading by own crown and shading by other trees (stand) separately and updates
// Qin_stand in TreeSegment.


//This functor EvaluateRadiationForCfTreeSegment evaluates shading
//caused by all other segments on this conifer segment. This functor
//uses functor ShadingEffectOfCfTreeSegment<TS,BUD> to go through all
//segments to check the shading.

//If the attributes voxel_space and border_forest are set (then
//voxel_space != NULL), the attenuation of the bean in the voxel_space
//and border_forest is taken into consideration.

template <class TS, class BUD, class TREE>
class EvaluateRadiationForCfTreeSegment_2 {
public:
    EvaluateRadiationForCfTreeSegment_2(const ParametricCurve& k, vector<TREE*>& vt):
        K(k), vtree(vt) {only_self = false;}
    EvaluateRadiationForCfTreeSegment_2(const ParametricCurve& k, vector<TREE*>& vt, bool o_s):
        K(k), vtree(vt), only_self(o_s) {}

    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;
private:
    const ParametricCurve& K;
    vector<TREE*>& vtree;
    bool only_self;
};


//===============================================================================
// VoxelSpace calculation

template <class TS, class BUD>
class EvaluateRadiationForCfTreeSegment_3 {
public:
    EvaluateRadiationForCfTreeSegment_3(const ParametricCurve& k) : K(k),
        evaluate_border(false) {}
    EvaluateRadiationForCfTreeSegment_3(const ParametricCurve& k,
                                        VoxelSpace* vs,BorderForest* bf, bool border,
                                        LGMdouble a, LGMdouble b, bool sd,LGMdouble kbc,
                                        bool d_e, bool wd, bool cs, LGMdouble st,bool calculateDirectionalStar):
        K(k),voxel_space(vs), border_forest(bf), evaluate_border(border),
        par_a(a), par_b(b), dump_self(sd), k_border_conifer(kbc), dir_effect(d_e),
        wood(wd), constant_star(st),correct_star(cs),calculateDirectionalStar(calculateDirectionalStar) {}

    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;

private:
    const ParametricCurve& K;
    VoxelSpace* voxel_space;
    BorderForest* border_forest;
    bool evaluate_border;
    LGMdouble par_a, par_b;
    bool dump_self;
    LGMdouble k_border_conifer;
    bool dir_effect;   //If effect of mean segment direction in boxes considered
    bool wood;         //If the needleless woody parts are included into calculation
    LGMdouble constant_star; //If this > 0, k = constant_star else k = star_mean
    bool correct_star;       //If star_eq -> star correction is done
    bool calculateDirectionalStar;



};

class AccumulateOpticalDepth{
 public:
 AccumulateOpticalDepth(LGMdouble side, LGMdouble a, LGMdouble b, Point loc, ParametricCurve kk,
			bool d_e, bool wd, bool cs, LGMdouble st,bool calculateDirectionalStar) :
  box_side_length(side), par_a(a), par_b(b), seg_loc(loc), K(kk), dir_effect(d_e), wood(wd), 
    constant_star(st),correct_star(cs), calculateDirectionalStar(calculateDirectionalStar)
  {box_volume = pow(box_side_length,3.0);
    int dummy = 0.0;
    //This function describes how much STAR deviates from mean STAR as a function of acute angle between
    //shoot and light beam directions. These values are in 
    //~/Riston-D/E/Hankkeet/LIGNUM/Erillishankkeet/Light/Radiation-article/mean-star-effect.dat". It has been
    //calculated in ~/Riston-D/E/Hankkeet/LIGNUM/Lignum-Forest/R-files/STAR.R
    ParametricCurve dir_effect_function(
"0 1.1629139 0.0785398163397448 1.14677595593143 0.15707963267949 1.13023014058558 0.235619449019234 1.11340055974152 "
"0.314159265358979 1.09637333348913 0.392699081698724 1.07916314315553 0.471238898038469 1.06142618194546 "
"0.549778714378214 1.04009823371317 0.628318530717959 1.01675794494258 0.706858347057703 0.991602507644704 "
"0.785398163397448 0.964551240973513 0.863937979737193 0.935664954902298 0.942477796076938 0.904595654266902 "
"1.02101761241668 0.871177309266682 1.09955742875643 0.83519042079002 1.17809724509617 0.795193529091556 "
"1.25663706143592 0.754853265395472 1.33517687777566 0.714203244987155 1.41371669411541 0.673172332712971 "
"1.49225651045515 0.631674241882798 1.5707963267949 0.5896551", dummy);
}

  double operator()(double o_d,VoxelMovement& vm){
    //    if((vm.af > R_EPSILON ||(wood && vm.wood_area > R_EPSILON)) && vm.n_segs_real > 0.0) {


    if(vm.af > R_EPSILON ||(wood && vm.wood_area > R_EPSILON)) {
      LGMdouble k = 0.14;   //= 0.14 is to avoid uninitialized warning message
      //**********Additional Code for the vectors of directional Star Values***************************************************
      vector<LGMdouble>  kdir(7);
      vector<LGMdouble>  angle(7);
      int count = 0;
      // Code added by KV this is the angle of inclination used for interpolation.
      for(double phi=0;phi<=PI_VALUE/2.0; phi+=PI_VALUE/12.0){
	angle[count] = phi;
	count+=1;
      }

      if(calculateDirectionalStar){
	kdir = vm.starDir;
      }
      else{
	if(constant_star > 0.0)
	  k = constant_star;
	else
	  k = vm.STAR_mean;
      }

      if(correct_star) {
	k = max(0.0,-0.014+1.056*k);
      }
      //****************************************************************************************************************
      //NOTE: here transformation STAR_eq --> STAR; documented in
      //~/Riston-D/E/LIGNUM/Light/summer-09-test/STAR-vs-STAR_eq.pdf

      //Effect of hit angle to the mean direction of shoots in the voxel box, documented in
      //~/Riston-D/E/LIGNUM/Light/Article/vs-STARmean-all.pdf and
      //~/Riston-D/E/LIGNUM/Light/Article/vs-STARmean-approximation.pdf
      PositionVector mean_dir = vm.mean_direction;
      LGMdouble mean_dir_length = mean_dir.length();
      LGMdouble effect = 1.0;
      if(dir_effect) {
	if(mean_dir_length > 0.0){
	  mean_dir.normalize();
	  LGMdouble inclination  =  PI_DIV_2 - acos(fabs(Dot(mean_dir,beam_dir)));
	  //	  effect =  K(inclination)/K(0.7);
	  effect = dir_effect_function(inclination);
	}
      }
      //this scales the effect depending on how parallel the segments are
      /*       if(mean_dir_length > 0.0) { */
      /* 	mean_dir.normalize(); */
      /* 	LGMdouble inclination  =  PI_DIV_2 - acos(fabs(Dot(mean_dir,beam_dir))); */
      /* 	//	LGMdouble effect  = 1.13 - 0.24*pow(inclination,2.0); */
      /* 	//	LGMdouble effect = 1.2-0.3*inclination*(inclination+0.3); */
      /* 	LGMdouble u =  mean_dir_length/vm.n_segs_real; */
      /* 	effect = 1.0 - u + u * K(inclination)/K(0.7); */
      /* 	cout << " u " << u << endl; */
      /*       } */

      //Code to evaluate k if calculateDirectionalStar is true then first search the lower bound then obtain the index by sebtracting the
      //beginning. Then add one to the lowest (or rather closest lower bound) index to get the upper bound index between which we need to
      //interpolate the new value.
      //Corrected by Risto 10.10.2014: sin(angle) and angle were compared. Now the calculation are made with angles.
      //Changes also in finding indexes
      if(calculateDirectionalStar){
	LGMdouble dir_angle = asin(beam_dir.getZ());
	std::vector<LGMdouble>::iterator up; //Up points to first element of angle which compares greater than dir_angle.
	up   = std::upper_bound(angle.begin(),angle.end(),dir_angle);
	if((up - angle.begin()) < 1) {
	  k = kdir[0];
	}
	else if((up - angle.begin()) > 6) {
	  k = kdir[6];
	}
	else {
	  int upperIndex = up - angle.begin();
	  int lowerIndex = upperIndex - 1;
	  k = kdir[lowerIndex] +(kdir[upperIndex]-kdir[lowerIndex])*
	    ((dir_angle-angle[lowerIndex])/(angle[upperIndex]-angle[lowerIndex]));
	  o_d += effect * k * vm.af * vm.l / box_volume;
	}
      }
      else{
	o_d += effect * k * vm.af * vm.l / box_volume;
      }
      if(wood) {
	//Mean projection area of surface of a circular cylinder (excluding end disks)
	// is 1/4 of its area
	o_d += 0.25 * vm.wood_area * vm.l / box_volume;
      }
      cout.precision(15);
    }
    return o_d;
  }

  //This is to take care of the direction of the beam of radiation
  PositionVector beam_dir;
 private:
  LGMdouble box_side_length;
  LGMdouble box_volume;
  LGMdouble par_a, par_b;
  Point seg_loc;   //location of segment
  ParametricCurve K;
  bool dir_effect;                  // If direction effect of segments in box considered
  bool wood;                       //  If woody parts are considered
  LGMdouble constant_star;        //   If this > 0, k = constant_star else k = star_mean
  bool correct_star;             //    If star_eq -> star correction is done
  bool calculateDirectionalStar;//     Directional Star values are required this is used
  ParametricCurve dir_effect_function;
};

//========================================================================================
//This functor is a modification of ShadingEffectOfCfTreeSegment_1<TS,BUD> and calculates transmission
//through a conifer segment to a point in space from different directions.
//The point and the directions are data members of this functor.
//The results are stored (and returned) in vector as optical depth S.

template <class TS,class BUD>
class ShadingEffectOfCfTreeSegmentToPoint {
public:
 ShadingEffectOfCfTreeSegmentToPoint(Point& p, vector<vector<LGMdouble> >& dir,
				     ParametricCurve& K_in, const vector<int>& ellipsoid_h,
                                     const bool ellipsoid_c, vector<LGMdouble>& S_in)
   : p0(p), directions(dir), K(K_in),ellipsoid_hits(ellipsoid_h),
    ellipsoid_calculation(ellipsoid_c), S(S_in)
  {
    number_of_directions = (int)directions.size();
    for(int i = 0; i < (int)directions.size(); i++)
      {S[i] = 0.0;}
  }
    TreeCompartment<TS,BUD>*  operator()(TreeCompartment<TS,BUD>* tc)const;
    vector<double>& getS(){return S;}
private:
    Point p0;
    vector<vector<LGMdouble> >& directions;
    const ParametricCurve& K;
    const vector<int>& ellipsoid_hits;
    const bool ellipsoid_calculation;      //
    vector<double>& S;           //Optical depth
    int number_of_directions;
};

void ellipsoid_interception(const Point& p00, const vector<vector<LGMdouble> >& directions,
			    const vector<vector<LGMdouble> >& ellipsoids,vector<int>& ellipsoid_hits);


#undef HIT_THE_FOLIAGE
#undef NO_HIT
#undef HIT_THE_WOOD

#include <CalculateLightI.h>

#endif
