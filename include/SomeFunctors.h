#ifndef SOME_FUNCTORS_H
#define SOME_FUNCTORS_H
#include <sstream>
#include <Lignum.h>
#include <ScotsPine.h>
#include <VoxelSpace.h>


///Propagate up the Qin to newly created segments and buds.
///Also set the relative light LGAip for the bud
class ForwardScotsPineQin{
public:
  double operator()(double& qin, TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      double ts_qin = GetValue(*ts,LGAQin);
      if (GetValue(*ts,LGAage) == 0.0){//pass qin to this segment
	SetValue(*ts,LGAQin,qin);
      }
      else{ //computed
	qin = ts_qin;
      }
    }
    else if (ScotsPineBud* b = dynamic_cast<ScotsPineBud*>(tc)){
      SetValue(*b,LGAQin,qin);
      Tree<ScotsPineSegment,ScotsPineBud>& t = 
	dynamic_cast<Tree<ScotsPineSegment,ScotsPineBud>&>(GetTree(*b));
      double max_light = GetValue(t,TreeQinMax);
      double ip = qin/max_light;
      SetValue(*b,LGAip,ip);
    }
    return qin;
  }
};

///Sets the variable apical, that creates the length difference between new segment
///on the basis of whether its apical (apical = 1) or lateral (apical = function of
///Qin of mother, that has been propagated to Qin of the new segment (by ForwardScotsPineQin())
class SetScotsPineSegmentApical{
 public:
  LGMdouble& 
    operator()(LGMdouble& go, TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const {

    if(ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if(GetValue(*ts,LGAage) == 0.0){
        if(GetValue(*ts, LGAomega) > go) {     //is lateral new segment
          ScotsPineTree& t = dynamic_cast<ScotsPineTree&>(GetTree(*ts));
          const ParametricCurve& fad = GetFunction(t,SPFAD);   //Function fapical.fun
          //Set the apical dominance as the function of relative light
	  double ap = fad(GetValue(*ts,LGAQin)/GetValue(t,TreeQinMax));
	  
	  // extern bool is_by_branches;
	  // if(is_by_branches) {
	  //     extern int ran3_seed;
	  //     ap *= (1.0 + 0.1*(ran3(&ran3_seed)-0.5));
	  //     if(ap > 1.0) ap = 1.0;
	  // }
          ts->setApical(ap);
        }
        else {
          ts->setApical(1.0);        //apical new segment
        }
        go = GetValue(*ts, LGAomega);  //This is in fact not necessary since new segment is last in axis      
      } else {
        go = GetValue(*ts, LGAomega);  //Old segment sends its Gravelius order
      }
    }
    return go;
  }
};



///Specific leaf  area for Scots pine  is a function  of relative light
///and set once for the newly created segments after the calculation of
///the light climate and the use of ForwardScotsPineQin. 
class SetScotsPineSegmentSf{
public:
  void operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0.0)
	SetValue(*ts,LGAsf);
    }
  }
};

///Set the needle angle once when  the segment is created as a function
///of light
class SetScotsPineSegmentNeedleAngle{
public:
  void operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0.0){
	const ScotsPineTree& t = dynamic_cast<ScotsPineTree&>(GetTree(*ts));
	const ParametricCurve& fna = GetFunction(t,SPFNA);
	double nl = GetValue(t,LGPnl);
	double B = GetValue(t,TreeQinMax);
	double qin = GetValue(*ts,LGAQin);
	//Update foliage height and radius to foliage limit
	SetValue(*ts,LGAHf,nl*sin(fna(qin/B)));
	SetValue(*ts,LGARf,GetValue(*ts,LGAR)+GetValue(*ts,LGAHf));
      }
    }
  }
};

/// \brief Collect the *sums* of main axis branch data
///
/// Collect the sums of branch diameters, branch lengths
/// and the number of branches to calculate branch means.
/// \sa Accumulate Branchmeans
class summing {
 public:
  summing():d2(0.0),d2l(0.0),lsum(0.0),n_br(0){}
  double d2;///< The diameter of branch (the first segment in Axis)
  double d2l;///< Defined as d2*L, L=branch length
  double lsum;///< Branch length as the sum of its segments
  int n_br; ///< Number of branches 
};
  
///\brief Calculate data for mean branch
class Branchmeans{
 public:
  /// The operator to be used  with Accumulate and summing
  /// \param id The identity element for Accumulate
  /// \param ts The tree segment for Accumulate
  /// \remark Implemented in branchfunctor.cc
  summing& operator ()(summing& id, TreeCompartment<ScotsPineSegment,
		       ScotsPineBud>* ts)const;
};

///\brief Set  STAR  mean based  on  segment age.
///
/// See  Smolander  et al  1994 TreePhys. Simply   initialize   the  functor
/// with   appropriate ParametricCurve and use ForEach.
template <class TS, class BUD>
class SetStarMean{
public:
  SetStarMean(const ParametricCurve f):fstar_mean(f){}
  TreeCompartment<ScotsPineSegment,ScotsPineBud>* operator() (TreeCompartment<TS,BUD>* tc) const{
    if (TS* ts = dynamic_cast<TS*>(tc)){
      LGMdouble starm = fstar_mean(GetValue(*ts,LGAage));
      SetValue(*ts,LGAstarm,starm);
    }
    return tc;
  }
private:
const ParametricCurve fstar_mean;
};





///A quick hack to try to  compute starSum and then star mean for Scots
///pine.
///
///This  is esseantially copy-paste from VoxelSpaceI.h  but I just
///take the  LGAstarm precomputed for  the segment to be  assigned into
///the voxel.   Quite a bit work for  small change.  What I  want to do
///(quickly) is to  compute the star mean for  the segment in different
///ways and add it to starSum. Then in the VoxelSpace::updateValues one
///can define starmean = starSum/Nsegments for each box.
inline
void DumpScotsPineSegment(VoxelBox &b, ScotsPineSegment& ts,double num_parts){	

  LGMdouble fmass = GetValue(ts, LGAWf) / num_parts;
  
  if (fmass > R_EPSILON){
    //LGAAf = LGPsf*LGAWf, and  LGPsf is now function of light
    LGMdouble farea = GetValue(ts,LGAAf)/num_parts;	
    b.addNeedleArea(farea);
    b.addNeedleMass(fmass);

    //Here  weight the  star mean  with  foliage area,  then in  the
    //VoxelBox::updateValues() remember to have: starm=starSum/AfTot
    b.addStarSum(GetValue(ts,LGAstarm)*farea);
    b.addWeight(farea);
    b.increaseNumberOfSegments();
  }
}

///Insert segment to  voxel. This is to compare  light methods 2
///and 4 with method 3.
inline 
void InsertSegmentToVoxel(VoxelBox &b, ScotsPineSegment& ts,double num_parts)
{
  LGMdouble fmass = GetValue(ts, LGAWf) / num_parts;
  
  if (fmass > R_EPSILON){
    //LGAAf = LGPsf*LGAWf, and  LGPsf is now function of light
    LGMdouble farea = GetValue(ts,LGAAf) / num_parts;	
    b.addNeedleArea(farea);
    b.addNeedleMass(fmass);

    //Here  weight the  star mean  with  foliage area,  then in  the
    //VoxelBox::updateValues() remember to have: starm=starSum/AfTot
    b.addStarSum(GetValue(ts,LGAstarm)*farea);
    b.addWeight(farea);
    
    //Qin and Qabs are interesting
    b.addQin(GetValue(ts,LGAQin)/num_parts);
    b.addQabs(GetValue(ts,LGAQabs)/num_parts);
    b.increaseNumberOfSegments();
  }
}

///This functor inserts  segments to voxels to compare  light methods 2
///and 4 with 3. Recall to  calculate mean qabs and mean qin by calling
///VoxelSpace::calculateMeanQabsQin().
template <class TS, class BUD>
class InsertSegmentsToVoxels{
public:
  InsertSegmentsToVoxels(VoxelSpace &s, int parts):space(s),num_parts(parts){}
  TreeCompartment<TS,BUD>* 
  operator ()(TreeCompartment<TS,BUD>* tc)const{
    if (TS* cfts = dynamic_cast<TS*>(tc)){
      if (GetValue(*cfts,LGAWf) > R_EPSILON){
	Point p = GetPoint(*cfts);
	PositionVector pv = GetDirection(*cfts);
	LGMdouble length = GetValue(*cfts, LGAL);
	//if the user wants 1 part (whole segment), the loop is executed
	//once  and the midpoint  of the  segment is  used; if  the user
	//wants 2 parts,  the loop is executed twice  and the points 1/3
	//and 2/3 of the segment length are used, and so on
	for (int i=1; i<(num_parts+1.0); i++){
	  Point p1 = p + (Point)((i/(num_parts+1.0))*length*pv);
          try{
            VoxelBox& box = space.getVoxelBox(p1);
            InsertSegmentToVoxel(box,*cfts, num_parts);
          }
          catch (OutOfVoxelSpaceException e){
            ostringstream error;
            Point b = e.getBox();
            Point p = e.getPoint();
            error << "DumpScotPineTree Functor Out of voxel space " << endl 
                  << "Box " << b.getX() << " " << b.getY() << " " << b.getZ() << " "
                  << "Point " << p.getX() << " " << p.getY() << " " << p.getZ() << endl;
            LGMError(error.str());
          } 
	}
	//Voxel space book keeping.  Descriptive data for Leaf Area Index:
	//add  needle  area.  Recall  after  last  insertion  to set  LAI:
	//Aftot/Avs (foliage area divided by voxel space area)
	ForestDescriptor& forest_descriptor = space.getForestDescriptor();
	SetValue(forest_descriptor,LGAAfc,
		 GetValue(forest_descriptor,LGAAfc)+GetValue(*cfts,LGAAf));
	//Crown limit update
	SetValue(forest_descriptor,LGAcbase,GetPoint(*cfts).getZ());
	//Crown heigth update
	SetValue(forest_descriptor,LGAH,GetPoint(*cfts,1.0).getZ());
      }
    } 
    return tc;
  }
private:
  VoxelSpace &space;
  double num_parts;//number of parts the segment is divided into (usually
  //1 but one may want to experiment)
};

template <class TS, class BUD>
class  DumpScotsPineTreeFunctor{
public:
  DumpScotsPineTreeFunctor(VoxelSpace &s, int parts):space(s),num_parts(parts){}  
  TreeCompartment<TS,BUD>* 
  operator ()(TreeCompartment<TS,BUD>* tc)const{
    if (TS* cfts = dynamic_cast<TS*>(tc)){
      if (GetValue(*cfts,LGAWf) > R_EPSILON){
	Point p = GetPoint(*cfts);
	PositionVector pv = GetDirection(*cfts);
	LGMdouble length = GetValue(*cfts, LGAL);
	//if the user wants 1 part (whole segment), the loop is executed
	//once  and the midpoint  of the  segment is  used; if  the user
	//wants 2 parts,  the loop is executed twice  and the points 1/3
	//and 2/3 of the segment length are used, and so on
	for (int i=1; i<(num_parts+1.0); i++){
	  Point p1 = p + (Point)((i/(num_parts+1.0))*length*pv);
          try{
            VoxelBox& box = space.getVoxelBox(p1);
            DumpScotsPineSegment(box,*cfts, num_parts);
          }
          catch (OutOfVoxelSpaceException e){
            ostringstream error;
            Point b = e.getBox();
            Point p = e.getPoint();
            error << "DumpScotPineTree Functor Out of voxel space " << endl 
                  << "Box " << b.getX() << " " << b.getY() << " " << b.getZ() << " "
                  << "Point " << p.getX() << " " << p.getY() << " " << p.getZ() << endl;
            LGMError(error.str());
          } 
	}
	//Voxel space book keeping.  Descriptive data for Leaf Area Index:
	//add  needle  area.  Recall  after  last  insertion  to set  LAI:
	//Aftot/Avs (foliage area divided by voxel space area)
	ForestDescriptor& forest_descriptor = space.getForestDescriptor();
	SetValue(forest_descriptor,LGAAfc,
		 GetValue(forest_descriptor,LGAAfc)+GetValue(*cfts,LGAAf));
	//Crown limit update
	SetValue(forest_descriptor,LGAcbase,GetPoint(*cfts).getZ());
	//Crown heigth update
	SetValue(forest_descriptor,LGAH,GetPoint(*cfts,1.0).getZ());
      }
    } 
    return tc;
  }
private:
  VoxelSpace &space;
  double num_parts;//number of parts the segment is divided into (usually
  //1 but one may want to experiment)
};

///Call this in the growth loop instead of DumpCfTree, num_parts is the
///number of parts  the segment is divided into (usually  1 but one may
///want to experiment)
inline void DumpScotsPineTree(VoxelSpace &s, 
			      Tree<ScotsPineSegment,ScotsPineBud> &tree,int num_parts){
  DumpScotsPineTreeFunctor<ScotsPineSegment,ScotsPineBud> f(s,num_parts);
  ForEach(tree, f);
}

inline
void SetScotsPineSegmentQabs(VoxelBox &b, ScotsPineSegment& ts, double num_parts)
{
  //Qin is needed for LGPsf
  SetValue(ts, LGAQin, GetValue(ts,LGAQin)+(b.getQin()/num_parts)); 
  //LGPAAf is LGPsf*LGAWf and LGPsf is now a function of light
  LGMdouble farea = GetValue(ts, LGAAf) / num_parts;
  LGMdouble qabs = 0.0;
  //Qabs computetd based on Qin, mean star and foliage area.
  qabs = b.getQin()*GetValue(ts,LGAstarm)*farea;
  if (qabs > b.getQin())
    cout << "QABS > QIN " <<  qabs << " " << b.getQin() << " "  << GetValue(ts,LGAstarm) << " "  << farea<< endl;
  //LGAAf, farea,  already takes account the num_parts
  SetValue(ts, LGAQabs,GetValue(ts,LGAQabs)+qabs);  
}

template <class TS, class BUD>
class SetScotsPineTreeQabsFunctor{
public:
  SetScotsPineTreeQabsFunctor(int parts):space(NULL),num_parts(parts){}
  TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;
  VoxelSpace* space;
private:
  double num_parts;
};

template <class TS, class BUD>
TreeCompartment<TS,BUD>* SetScotsPineTreeQabsFunctor<TS,BUD>::operator ()(TreeCompartment<TS,BUD>* tc)const
{
  if (TS* cfts = dynamic_cast<TS*>(tc))
    {
      LGMdouble fmass = GetValue(*cfts, LGAWf);
      //Reset Qabs Qin  here
      SetValue(*cfts, LGAQabs, 0.0);
      SetValue(*cfts, LGAQin, 0.0);
      if (fmass > R_EPSILON)
	{
	  Point p = GetPoint(*cfts);
	  PositionVector pv = GetDirection(*cfts);
	  LGMdouble length = GetValue(*cfts, LGAL);
	  //if the user wants 1 part (whole segment), the loop is executed
	  //once  and the midpoint  of the  segment is  used; if  the user
	  //wants 2 parts,  the loop is executed twice  and the points 1/3
	  //and 2/3 of the segment length are used, and so on
	  for (int i=1; i<(num_parts+1.0); i++)
	    {
	      Point p1 = p + (Point)(length * (i/(num_parts+1.0)) * pv);
	      VoxelBox box = space->getVoxelBox(p1);
	      SetScotsPineSegmentQabs(box, *cfts, num_parts);
	    }		
	}
    }
  return tc;    
}

template <class TS,class BUD>
void SetScotsPineTreeQabs(VoxelSpace &s, Tree<TS, BUD> &tree, int num_parts)
{
  SetScotsPineTreeQabsFunctor<TS,BUD> f(num_parts);
  f.space = &s;
  ForEach(tree, f);
}

    class SetSapwoodDemandAtJunction {
  public:
      double& operator()(double& oomega, TreeCompartment<ScotsPineSegment,
			 ScotsPineBud>* tc)const
      {
        if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment *>(tc))
          {
            double myomega = GetValue(*ts,LGAomega);
            if(myomega - R_EPSILON  > oomega) {
              const ParametricCurve& pc =
 		GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPSD);
              SetValue(*ts,SPAAsDown,pc(oomega));
              oomega = myomega;
            }
            else
              SetValue(*ts,SPAAsDown,1.0);
          }
        return oomega;
      }
  };

class BranchInformation {
public:
  BranchInformation(const string& filename): fname(filename)
  {
    ofstream f(fname.c_str() , ofstream::trunc);
    f << "Height Dist_top Diameter HwDiam Length Foliage_mass  Sapwood_mass Heartwood_mass Wood_mass"
         " Photosynthesis Respiration rel_pos" << endl;
    f.close();
  }

  TreeCompartment<ScotsPineSegment,ScotsPineBud>*
		operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
  {
    if (Axis<ScotsPineSegment,ScotsPineBud>*ax =
	dynamic_cast<Axis<ScotsPineSegment,ScotsPineBud>*>(tc)){
     list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>& cmpls =
	GetTreeCompartmentList(*ax);
      if(cmpls.size() > 2) {
	TreeSegment<ScotsPineSegment,ScotsPineBud>* fs = GetFirstTreeSegment(*ax);
	if((GetValue(*fs, LGAomega) > 1.0) && (GetValue(*fs, LGAomega) < 3.0)) {
	  Tree<ScotsPineSegment,ScotsPineBud>& tt = GetTree(*fs);
	  LGMdouble h0 = GetValue(tt,LGAH);
	  LGMdouble hh = GetPoint(*fs).getZ();
	  LGMdouble d = 200.0*GetValue(*fs, LGAR);
	  LGMdouble dh = 200.0*GetValue(*fs, LGARh);
	  LGMdouble Wf = 2.0*GetBranchFoliage(*ax);
	  LGMdouble Ws = 2.0*GetBranchSapwoodMass(*ax);
	  LGMdouble Wh = 2.0*GetBranchHeartwoodMass(*ax);
	  LGMdouble Ww = 2.0*GetBranchWoodMass(*ax);
	  LGMdouble P = 2.0*GetBranchPhotosynthesis(*ax);
	  LGMdouble M = 2.0*GetBranchRespiration(*ax);

	  //Relative position of the branch (0 = top, 1 = crown base)
	  Tree<ScotsPineSegment,ScotsPineBud>& tree = GetTree(*fs);
	  DCLData dcl;
	  AccumulateDown(tree,dcl,AddBranchWf(),DiameterCrownBase<ScotsPineSegment,ScotsPineBud>());
	  LGMdouble Hc = dcl.HCrownBase();
	  LGMdouble rel_pos = (h0 - hh) /(h0 - Hc);
	  rel_pos = max<double>(0.0,min<double>(1.0,rel_pos));

	  ofstream f(fname.c_str() , ofstream::app);
	  f << hh << " " << h0 - hh << " " << d << " " << dh << " " << GetValue(*ax,LGAL) << " "
	    << Wf << " " << Ws << " " << Wh << " " << Ww << " " << P << " " << M << " " << rel_pos << endl;
	  f.close();
	}
      }

    }
    return tc;
  }

 private:
  string fname;
};



///Evaluate the surface area of the woody part (excluding end discs) of
///new segments (that have foliage)
class SurfaceAreaOfNewSegments{
public:
  double operator()(double& A, TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) > R_EPSILON && GetValue(*ts,LGAWf) > R_EPSILON) {
	double L = GetValue(*ts, LGAL);
	double R = GetValue(*ts, LGAR);
	A += 2.0 * PI_VALUE * R * L;
      }
    }
    return A;
}

};

///Vertical distribution of  fip(ip) in segments The vector  v has Tree
///age  positions  representing  the  vertical  division  of  the  tree
///into intervals  of mean annual growth).
///Usage is simply:<br>
///   pair<vector<pair<double,int> >,double> p;<br>
///   Accumulate(tree,v,CollectVerticalDistributionOfFip());<br>
///where the  p.first is the vector  for cumulative fips  and number of
///their observations and p.second the height interval When printing to
///the  file  remember  to  divide  cumulative "fips"  with  number  of
///observations.
class CollectVerticalDistributionOfFip{
public:
  pair<vector<pair<double,int> >,double>& 
  operator()(pair<vector<pair<double,int> >,double>& p,
	     TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAWf) > R_EPSILON){
	const ParametricCurve& fip = GetFunction(GetTree(*ts),LGMIP);
	double B = GetValue(GetTree(*ts),TreeQinMax);
	double qin = GetValue(*ts,LGAQin);
	double ip = qin/B;
	double result = fip(ip);
	double height = GetPoint(*ts).getZ();
	//Find the height category for the segment (p.second is the interval)
	double category = height/p.second;
	//Drop the decimal part, the integer part is the index for the
	//position in the vector v.
	unsigned int index = static_cast<unsigned int>(category);
	//cout << category << " " << index << " "  << height << " "  << p.second << endl;
	if (index < p.first.size()){
	  p.first[index].first += result;//cumulative f(ip) 
	  p.first[index].second += 1; //number of  observations (segments with
	                              //foliage)
	}
	else{
	  cout << "CollectVerticalDistributionOfFip: Segment height exceeds tree height" <<endl;
	  cout << "Ignoring " << height << " > "  << GetValue(GetTree(*ts),LGAH);
	}
      }
    }
    return p;
  }
};

class CollectSegmentDiameters{
public:
  list<double>& operator()(list<double>& ls,TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      double d = 2.0* GetValue(*ts,LGAR);
      ls.push_back(d);
    }
    return ls;
  }
};

class SegmentProductionBalance {
 public:
  SegmentProductionBalance(const string& filename): fname(filename)
    {
      ofstream f(fname.c_str() , ofstream::trunc);
      f << "X Y Z Dist_top Dist_stem Vigor Diam Length   Age    Wf   Qin   Qabs   Assim"
	   "   Respiration" << endl;
      f.close();
    }

    TreeCompartment<ScotsPineSegment,ScotsPineBud>*
      operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
      {
	if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	  LGMdouble Wf = GetValue(*ts,LGAWf);
	  if(Wf > 0.0) {
	    LGMdouble P = GetValue(*ts,LGAP);
	    LGMdouble M = GetValue(*ts,LGAM);
	    LGMdouble Qin = GetValue(*ts,LGAQin);
	    LGMdouble Qabs = GetValue(*ts,LGAQabs);
	    LGMdouble D = 200.0*GetValue(*ts,LGAR);
	    LGMdouble L = 100.0*GetValue(*ts,LGAL);
	    LGMdouble Age = GetValue(*ts, LGAage);
	    LGMdouble Vi = GetValue(*ts, LGAvi);
	    Tree<ScotsPineSegment,ScotsPineBud>& t = GetTree(*ts);
	    LGMdouble H = GetValue(t, LGAH);
	    Point pp = GetPoint(*ts);
	    LGMdouble ver_dist = H - pp.getZ();
	    Axis<ScotsPineSegment,ScotsPineBud>& ax = GetAxis(t);
	    TreeCompartment<ScotsPineSegment,ScotsPineBud>* fc = GetFirstTreeCompartment(ax);
	    LGMdouble hor_dist = sqrt(pow(pp.getX()-GetPoint(*fc).getX(),2.0) +
				      pow(pp.getY()-GetPoint(*fc).getY(),2.0));
	
	    ofstream f(fname.c_str() , ofstream::app);
	    PositionVector mp = PositionVector(pp)+0.5*GetValue(*ts,LGAL)*GetDirection(*ts);
	    f << pp.getZ() << mp.getX() << " " << mp.getY() << " "
	      << mp.getZ() << " " << ver_dist << " " << hor_dist << " " << Vi << " " << D << " " << L << " "
	      << Age << " " << 2000.0*Wf << " " << Qin << " " << Qabs << " " << 2000.0*P << " " << 2000.0*M
	      << endl;
	    f.close();
	  }
	}
	return tc;
      }

 private:
    string fname;
};


///Evaluate the surface area of the woody part (excluding end discs) of
///new segments (that have foliage)
class SetRh{
 public:
 SetRh(const LGMdouble pc) : per_cent(pc) {}
  TreeCompartment<ScotsPineSegment,ScotsPineBud>* operator()
    (TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
    {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	LGMdouble r = GetValue(*ts,LGAR);
	SetValue(*ts, LGARh, per_cent*r);
      }
      return tc;
    }
 private:
  LGMdouble per_cent;
};


class Sum2ndOrderBranchesSegments {
 public:
  void operator()(Axis<ScotsPineSegment,ScotsPineBud>& axis, vector<int>& values)const
  {
    list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>& tc_ls =
      GetTreeCompartmentList(axis);
    list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>::iterator i_tcl;
    for (i_tcl = tc_ls.begin(); i_tcl != tc_ls.end();i_tcl++){
      if(ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(*i_tcl)) {
        values[1] += 1;
        if(GetValue(*ts, LGAWf) > R_EPSILON) {
          values[2] += 1;
        }
      }
      else if(BranchingPoint<ScotsPineSegment,ScotsPineBud>* bp =
              dynamic_cast<BranchingPoint<ScotsPineSegment,ScotsPineBud>*>(*i_tcl)){
        list<Axis<ScotsPineSegment,ScotsPineBud>*>& axis_ls = GetAxisList(*bp);
        list<Axis<ScotsPineSegment,ScotsPineBud>*>::iterator i_axl;
        for (i_axl = axis_ls.begin(); i_axl != axis_ls.end(); i_axl++) {
	  if(GetTreeCompartmentList(**i_axl).size() > 2) {
	    Bud<ScotsPineSegment, ScotsPineBud>* tb = GetTerminatingBud(**i_axl);
            if(GetValue(*tb,LGAomega) == 3.0) {
              values[0] += 1;
            }
	    vector<int> v1(3,0);
	    (*this)(**i_axl, v1);
	    values[0] += v1[0];
	    values[1] += v1[1];
	    values[2] += v1[2];
	  }
        } //for (i_axl = axis_ls.begin(); i_axl ...
      }
    }
    return;
  }
};


///UnDumpScotsPineSegment  is   the  inverse  to  DumpScotsPineSegment,
///instead of adding foliage remove it
inline void UnDumpScotsPineSegment(VoxelBox &b, ScotsPineSegment& ts,double num_parts,
				   bool wood_voxel){	

  LGMdouble fmass = GetValue(ts, LGAWf) / (double)num_parts;  
  LGMdouble farea = GetValue(ts,LGAAf) / (double)num_parts;	
  b.subtractNeedleArea(farea);
  b.subtractNeedleMass(fmass);

  //Here  weight the  star mean  with  foliage area,  then in  the
  //VoxelBox::updateValues() remember to have: starm=starSum/AfTot
  b.subtractStarSum(GetValue(ts,LGAstarm)*farea);
  b.subtractWeight(farea);

  b.decreaseNumberOfSegments();

  b.addNumberOfSegmentsReal(-1.0/(double)num_parts); 

  if(wood_voxel) {
      LGMdouble r = GetValue(ts, LGAR);
      LGMdouble length = GetValue(ts, LGAL) / (double)num_parts;
      LGMdouble mass = GetValue(ts, LGAWood) / (double)num_parts;
      LGMdouble area = 2.0*PI_VALUE*r*length;
      b.subtractWoodMass(mass);
      b.subtractWoodArea(area);
   }
 }



template <class TS, class BUD>
  class  UnDumpScotsPineTreeFunctor{
 public:
  UnDumpScotsPineTreeFunctor(VoxelSpace &s, int parts, bool wood_vox):space(s),num_parts(parts),
    wood_voxel(wood_vox){}  
    TreeCompartment<TS,BUD>* 
      operator ()(TreeCompartment<TS,BUD>* tc)const{
      if (TS* cfts = dynamic_cast<TS*>(tc)){
	bool foliage = false;
	if (GetValue(*cfts,LGAWf) > R_EPSILON)
	  foliage = true;
	if(foliage || wood_voxel) {
	  Point p = GetPoint(*cfts);
	  PositionVector pv = GetDirection(*cfts);
	  LGMdouble length = GetValue(*cfts, LGAL);
	  //if the user wants 1 part (whole segment), the loop is executed
	  //once  and the midpoint  of the  segment is  used; if  the user
	  //wants 2 parts,  the loop is executed twice  and the points 1/3
	  //and 2/3 of the segment length are used, and so on
	  for (int i=1; i<(num_parts+1.0); i++){
	    Point p1 = p + (Point)(((double)i/((double)num_parts+1.0))*length*pv);
	    try{
	      VoxelBox& box = space.getVoxelBox(p1);
	      UnDumpScotsPineSegment(box,*cfts, num_parts,wood_voxel);
	    }
	    catch (OutOfVoxelSpaceException e){
	      ostringstream error;
	      Point b = e.getBox();
	      Point p = e.getPoint();
	      error << "DumpScotPineTree Functor Out of voxel space " << endl 
		    << "Box " << b.getX() << " " << b.getY() << " " << b.getZ() << " "
		    << "Point " << p.getX() << " " << p.getY() << " " << p.getZ() << endl;
	      LGMError(error.str());
	    }
  
	  }  //for(int i = 1 ....
	} //if(foliage || ...
      } 
      return tc;
    }
 private:
    VoxelSpace &space;
    double num_parts;//number of parts the segment is divided into (usually
    //1 but one may want to experiment)
    bool wood_voxel; //If woody part is considered also
};


///UnDumpScotsPineTree is the "inverse" to DumpCfTree in VoxelSpaceI.h and VoxelVoxI.h,
///it removes the values of this tree from the voxels
inline void UnDumpScotsPineTree(VoxelSpace &s, 
				Tree<ScotsPineSegment,ScotsPineBud> &tree,int num_parts,
				bool wood_voxel){
  UnDumpScotsPineTreeFunctor<ScotsPineSegment,ScotsPineBud> f(s,num_parts,wood_voxel);
  ForEach(tree, f);
}

///Sort a vector containing pairs `p`, where p.first is the tree position in `vtree`
///and `p.second` the tree height, in ascending order based on tree height
class SortByTreeHeight{
public:
  bool operator()(pair<unsigned int,double>& p1, pair<unsigned int,double>& p2)
  {
    if (p1.second < p2.second){
      return true;
    }
    else{
      return false;
    }
  }
};

class CollectTreePositions{
public:
  unsigned int operator()(const pair<unsigned int,double>& p){
    return p.first;
  }
};
#endif
