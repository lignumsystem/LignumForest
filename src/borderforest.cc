#include <BorderForest.h>
#include <Nearby.h>
namespace Lignum {

  //Return the extinction caused by the border stand
  //Input: p0   start point of the light beam
  //       dir  direction of the light beam, |dir| == 1 (!!!)
  //Calculate the  point where  the light beam  exits the  voxel space
  //(there  must be  one). NearByShading  then returns  the extinction
  //coeffcient
  double BorderForest::getBorderForestExtinction(const Point& p0, const PositionVector& dir,
						 LGMdouble k_conifer)
  {
    //Start point of the light beam
    PositionVector d0(p0);
    //Normals to the faces of the voxel space
    //PositionVector n1(1,0,0);//normal of the front face of the voxel space
    //PositionVector n2(0,1,0);//normal of the left face of the voxel space
    //PositionVector n3(0,0,1);//normal of the bottom face of the voxel space
    //PositionVector n4(-1,0,0);//normal of the back face of the voxel space
    //PositionVector n5(0,-1,0);//normal of the right face of the voxel space
    //PositionVector n6(0,0,-1);//normal of the top face of the voxel space
    //Origo  of  the  voxel  space in  global  (segment)  coordinates,
    //i.e. the point on the front,  left and bottom faces of the voxel
    //space
    //Point p1(corner1);
    Point p1(corner_l.getX(),corner_l.getY(),0.0);   // Works only for stand and
                       // border forest standing on level  z = 0

    //opposite point  to origo  in global (segment)  coordinates, i.e.
    //the point on the back, right and top faces of the voxel space
    //Point p2(corner2);

    Point p2(corner_r.getX(), corner_r.getY(), H);

    //Calculate the  distances light  beam can travel  before crossing
    //the voxel space in x,y and  z directions. This is the problem of
    //deciding  if  a  ray  intersects  with  a  plane.   The  ray  is
    //represented as d0+t*dir, where d0  is the starting point and dir
    //is the direction  (unit vector) of the ray.  't' is the distance
    //to the plane. The plane is represented as Ax+By+Cz+D=0, where A,
    //B and C  is the normal to  the plane (unit vector) and  D is the
    //(shortest)  distance of  the plane  to  origo. At  the point  of
    //intersection   the    ray   satisfies   the    plane   equation:
    //A*(d0.x+t*dir.x)+B*(d0.y+t*dir.y)+C*(d0z+t*dir.z)+D=0  Solve the
    //equation for t:
    //t=-(A*d0.x+B*d0.y+C*d0.z+D)/(A*dir.x+B*dir.y+C*dir.z)  Note  the
    //sign of D; it is a positive number in Ax+By+Cz=D and negative in
    //Ax+By+Cz+D=0. Note also that the  normals are simple and we know
    //the D, so the equation for t simplifies quite a lot.
    double t1,t2,t3,t4,t5,t6;
    t1=t2=t3=t4=t5=t6=-1.0;//initialize to negative (i.e. no  intersection)
    if (fabs(dir.getX()) > R_EPSILON){
      //A=1,B=C=0,D=-corner1.X 
      t1 = -(d0.getX() + (-p1.getX()))/(dir.getX());//front face
      //A=1,B=C=0,D=-corner2.X
      t4 = -(d0.getX() + (-p2.getX()))/(dir.getX());//back face
    }
    if (fabs(dir.getY()) > R_EPSILON){
      t2 = -(d0.getY() + (-p1.getY()))/(dir.getY());//left face
      t5 = -(d0.getY() + (-p2.getY()))/(dir.getY());//right face
    }
    if  (fabs(dir.getZ()) > R_EPSILON){
      t3 = -(d0.getZ() + (-p1.getZ()))/(dir.getZ());//bottom face
      t6 = -(d0.getZ() + (-p2.getZ()))/(dir.getZ());//top face
    }
    vector<double> v(6,0.0);
    v[0] = t1; v[1] = t2; v[2] = t3; v[3] = t4; 
    v[4] = t5; v[5] = t6;
    //Sort in ascending order
    sort(v.begin(),v.end());
    //Take  the first nonnegative  t, i.e.  the shortest  distance the
    //beam can travel in the voxel space before crossing some wall
    vector<double>::iterator it = find_if(v.begin(),v.end(),
					  bind2nd(greater_equal<double>(),0.0));
    double tdist = R_HUGE;
    if (it == v.end()){
      cerr << "No Exit point from voxel space (All t < 0). Error!!!" << endl;
      cerr << "Start point " << d0 <<endl;
      cerr << "Beam direction " << dir <<endl;
      cerr << "t1,....,t6 " << flush;
      copy(v.begin(),v.end(),ostream_iterator<double>(cerr," "));
      cerr << endl;
    }
    else{
      tdist = *it;
    }
    //The exit point from the voxel space

    //Now: The exit point from the stand space to border forest
    PositionVector exit = d0+tdist*dir;
//     double tau = NearbyShading(Point(exit),dir,
// 			       GetValue(forest_descriptor,LGAH),
// 			       GetValue(forest_descriptor,LGAcbase),
// 			       GetValue(forest_descriptor,LGALAIc),
// 			       GetValue(forest_descriptor,LGALAIb));

    LGMdouble k_deciduous = 0.0;     //This time
    double tau = NearbyShading(Point(exit),dir,
			       H, Hcb, LAI, 0.0, k_conifer, k_deciduous);     //LAI of deciduous for the time being = 0
 
    return tau;
  }
    
 
}  //Namespace Lignum
