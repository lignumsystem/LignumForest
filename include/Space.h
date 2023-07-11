

#ifndef SPACE_H
#define SPACE_H
#include <sstream>
#include <Lignum.h>
#include <ScotsPine.h>
#include <VoxelSpace.h>
///\file Space.h
///\brief This file contains stuff about space colonization etc.

namespace LignumForest{
  ///\brief https://en.wikipedia.org/wiki/Spherical_sector
  ///
  ///In geometry, a spherical sector is a portion of a sphere defined by
  ///a conical boundary with apex at the center of the sphere. It can be described
  ///as the union of a spherical cap and the cone formed by the center of the sphere
  ///and the base of the cap.
  class SphericalSector {
  public:
    SphericalSector(const Point& iapex, const PositionVector& idirection,
		    const LGMdouble& ihalf_angle, const LGMdouble& iheight) :
      apex(iapex), direction(idirection), half_angle(ihalf_angle),
      height(iheight) {direction.normalize(); /*precaution*/}
    void getExtremePoints(const int& n_rim, list<Point>& cone_points) {
      //returns apex and n_rim points at a circle at distance r from
      //the apec so that the points on the circle "enclose" the spehrical
      //cap

      cone_points.push_back(apex);
      Point rim_base = apex + height * Point(direction);
      PositionVector end_v;    //vector perpendicular to direction
      if((direction || PositionVector(0.0,0.0,1.0)) > R_EPSILON){
	end_v = Cross(PositionVector(0.0,0.0,1.0),direction);
      }
      else {
	end_v = Cross(PositionVector(1.0,0.0,0.0),direction);
      }
      end_v.normalize();
      end_v = height * tan(half_angle) * end_v;
 
      double r_angle = 2.0*PI_VALUE/static_cast<LGMdouble>(n_rim);
      for(int i = 0; i < static_cast<LGMdouble>(n_rim); i++) {
	Point p_rim = rim_base + Point(end_v);
	cone_points.push_back(p_rim);
	end_v.rotate(PositionVector(0.0,0.0,0.0), direction, r_angle);
      }
    }
 
    bool in_SphericalSector(const Point& p);
  private:
    Point apex;
    PositionVector direction;
    LGMdouble half_angle;
    LGMdouble height;
  };


  class SetBudViewFunctor
  {
  public:
    SetBudViewFunctor(  VoxelSpace *ivs, const LGMdouble& ch, const  LGMdouble& cha, const int& rp) : 
      vs(ivs), cone_h(ch), cone_ha(cha), points_on_rim(rp) {
      cone_volume = 2.0*PI_VALUE*pow(cone_h,3.0)*(1.0-cos(cone_ha))/3.0;}

    TreeCompartment<ScotsPineSegment,ScotsPineBud>* 
    operator ()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const {
      if (ScotsPineBud* bud = dynamic_cast<ScotsPineBud*>(tc)){
	if(GetValue(*bud, LGAstate) == ALIVE) {

	  //First cone points into list: apex, and a number of points
	  //in the circular outer rim
	  Point p0 = GetPoint(*bud);
	  PositionVector b_dir = GetDirection(*bud);
	  SphericalSector ss(p0, b_dir, cone_ha, cone_h);
	  list<Point> cone_points;
	  ss.getExtremePoints(points_on_rim, cone_points);

	  //Now bounding box around the cone
	  BoundingBox cone_bbox;
	  list<Point>::iterator I;
	  for(I = cone_points.begin(); I != cone_points.end(); I++) {
	    cone_bbox.addPoint(*I);
	  }
	  vector<int> lli = vs->getBoxIndexes(cone_bbox.getMin());
	  //	cout << "ll " << lli[0] << " " << lli[1] << " " << lli[2] << endl;
	  vector<int> uri = vs->getBoxIndexes(cone_bbox.getMax());
	  //	cout << "ur " << uri[0] << " " << uri[1] << " " << uri[2] << endl;

	  //	cout << "Max I " << vs->getNoBoxX() << " " << vs->getNoBoxY() << " " << vs->getNoBoxZ()<< endl;
	  int min_ix = max(0,lli[0]);
	  int min_iy = max(0,lli[1]);
	  int min_iz = max(0,lli[2]);
	  int max_ix = min(vs->getNoBoxX(),uri[0]);
	  int max_iy = min(vs->getNoBoxY(),uri[1]);
	  int max_iz = min(vs->getNoBoxZ(),uri[2]);

	  /* 	cout << "min " << min_ix << " " << min_iy << " " << min_iz  << " " << endl; */
	  /* 	cout << "max " << max_ix << " " << max_iy << " " << max_iz  << " " << endl; */



	  LGMdouble af_sum = 0.0;
	  for(int i = min_ix; i < max_ix; i++) {
	    for(int j = min_iy; j < max_iy; j++) {
	      for(int k = min_iz; k < max_iz; k++) {
		//	      cout << i << " " << j << " " << k << " " << vs->voxboxes[i][j][k].getCenterPoint();
		if(ss.in_SphericalSector(vs->voxboxes[i][j][k].getCenterPoint())) {
		  af_sum += vs->voxboxes[i][j][k].getNeedleArea();
		}
	      }
	    }
	  }


	  LGMdouble view = 0.0;
	  view = af_sum/cone_volume;

	  bud->view = view;
	}
      }
      return tc;
    }

  private:
    VoxelSpace *vs;
    LGMdouble cone_h;    //cone height [m]
    LGMdouble cone_ha;   //cone half angle
    int points_on_rim;   //no points on spherical sector out rim for bounding box
    LGMdouble cone_volume;
	
  };
}//End namespace LignumForest
#endif
