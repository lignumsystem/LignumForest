#ifndef SCOTSPINE_H
#define SCOTSPINE_H
#include <Pine.h>
#include <VoxelSpace.h>

//SPFN =  Scots Pine Functions,
//SPFAF = Initial foliage m^2/kgC   
//SPFAD = Scots Pine  Apical Dominance,
//SPFGO = Scots Pine Gravelius Order
//SPFLR = Scots Pine length - radius relationship, R = f(relative_light)*R 
//SPFNA = Scots Pine Needle Angle
//SPFSF = Scots Pine Specific Leaf Area
//SPSD = Scots Pine Sapwood Down as a function of Gravelius order
//SPEBHF = Scots Pine Extended Borcher-Honda lambda value function
//SPBVF = Scots Pine Bud View function
enum SPFN {SPFAF,SPFAD, SPFGO,SPFLR,SPFNA,SPFSF, SPSD, SPWD, SPEBHF, SPBVF};

//Enumeration for SetValue, GetValue in ScotsPine
//Scots Pine Attribute Double SPAD,
//Sapwood down, Height at crown limit, start of heartwood build up 
enum SPAD {SPAAsDown,SPCrownRatio,SPHc,SPHwStart, SPrue};

//Scots Pine Parameter Double SPPD
//Extended Borchert-Honda (1 = true,  < 1 = false)
enum SPPD {SPis_EBH};

// 0  LGAplength  Path length  from the base of the  tree to a segment        
class ScotsPineBud;
class ScotsPineSegment;

extern ParametricCurve bud_view_f;

extern ParametricCurve adhoc;
extern bool is_adhoc;
extern VoxelSpace space_occupancy;

extern bool space0;
extern bool space1;
extern bool space2;
extern double space2_distance;

extern double L_age, L_H;
extern double global_hcb;

class ScotsPineTree: public Tree<ScotsPineSegment,ScotsPineBud>{
  //Return the function asked
  friend const ParametricCurve& GetFunction(const ScotsPineTree& t, SPFN name)
  {
    switch (name){
    case SPFAF:
      return t.af;//Intial foliage m^2/kgC as a function of light
    case SPFAD:
      return t.adf;//Apical dominance as a function of light
    case SPFGO:
      return t.gof;//Gravelius order effect on segment length 
    case SPFLR:
      return t.lr;//Length radius relationship as a function of light
    case SPFNA:
      return t.naf;//Needle angle as a function of light
    case SPFSF://Specific leaf area as a function of light
      return t.sf;
    case SPSD: //Scots Pine Sapwood Down as a function of Gravelius order
      return t.spsd;
    case SPWD:
      return t.wd;
    case SPEBHF: //Scots Pine Extended Borchert-Honda lambda value function
      return t.ebhf;
    case SPBVF: //Scots Pine Bud View function
      return t.bvf;
    default:
      LGMMessage("ScotsPineTree Unknown function");
      throw ParametricCurve();
    }
  }


  ///Two alternatives of SetFunction, the latter one comes from merging tree
  ///growth of the project CrownDensity
  friend void SetFunction(ScotsPineTree& t, ParametricCurve& f, SPFN name)
  {  
    if (name == SPFAF){
      t.af = f;
    }
    else if (name == SPFAD){
      t.adf = f;
    }
    else if (name == SPFGO){
      t.gof = f;
    }
    else if (name == SPFLR){
      t.lr = f;
    }
    else if (name == SPFNA){
      t.naf = f;
    }
    else if (name == SPFSF){
      t.sf = f;
    }
    else if (name == SPSD){
      t.spsd = f;
    }
    else if (name == SPWD){
      t.wd = f;
    }
    else if (name == SPEBHF){
      t.ebhf = f;
    }
    else{
      cerr << "SetFunction unknown function: " << name << endl;
      throw ParametricCurve();
    }
    return;
  }

  // friend void SetFunction(ScotsPineTree& t, SPFN name, ParametricCurve funct)
  // {
  //   cout << name << endl;
  //   switch (name){
  //   case SPEBHF: //Scots Pine Extended Borchert-Honda lambda value function
  //     {
  // 	(t.ebhf).install(funct.getFile());
  //     }
  //   default:
  //     LGMMessage("ScotsPineTree Unknown function name");
  //     throw ParametricCurve();
  //   }
  // }
  

  friend double GetValue(const  ScotsPineTree& t, SPAD name)
  {
    if (name == SPHwStart){
      return t.hw_start;
    }
    else if (name == SPHc){
      return t.hc;
    }
    else if (name == SPCrownRatio){//crown ratio
      return (GetValue(t,LGAH)-GetValue(t,SPHc))/GetValue(t,LGAH);
    }
    else{
      cerr << "GetValue(ScotsPineTree,name) unknown name " << name <<endl;
      return t.hw_start;
    }
  }

  friend double GetValue(ScotsPineTree& t, SPPD name) {
    if (name == SPis_EBH){ //If Palubicki et al 2009 Extended Borchert et Honda
      // resource distribution in use (0 = false, 1 = true)
      return t.is_EBH;
    }
    else{
      cerr << "GetValue(ScotsPineTree,SPPD name) unknown SPPD parameter name " << name <<endl;
      return 0.0;
    }
  }


  friend double SetValue(ScotsPineTree& t, SPAD name, double value)
  {
    double old_value = GetValue(t,name);
    if (name == SPHwStart){
      t.hw_start = value;
    }
    else if (name == SPHc){//height at crown limit
      t.hc = value;
    }
    else{
      cerr << "GetValue(ScotsPineTree,name) unknown SPAD name " << name <<endl;
    }
    return old_value;
  }

  friend double SetValue(ScotsPineTree& t, SPPD name, double value)
  {
    //    double old_value = GetValue(t,name);
    if (name == SPis_EBH){//if Extended Borchert et Honda distn in use
      t.is_EBH = value;        //(0 = false, 1 = true)
    }
    else{
      cerr << "GetValue(ScotsPineTree, SPPD name) unknown name " << name <<endl;
    }
    return 0.0;
  }

  
public:
  ScotsPineTree(const Point& p, const PositionVector& d,
		const string& sffun, const string& apicalfun, 
		const string& gofun, const string& sapwdownfun,
		const string& affun, const string& nafun, 
		const string& wdfun, const string& lrfun, const string& ebhfun,
		const string& bvffun)
    :Tree<ScotsPineSegment,ScotsPineBud>(p,d),sf(sffun),adf(apicalfun),
     gof(gofun),spsd(sapwdownfun),af(affun),naf(nafun),wd(wdfun),lr(lrfun),
     ebhf(ebhfun),bvf(bvffun),hw_start(0.0),hc(0.0),is_EBH(0.0){}
  LGMdouble getBranchAngle() {return branch_angle;}
  void setBranchAngle(LGMdouble ba) {branch_angle = ba;}

private:
  ParametricCurve sf;  ///Specific leaf area: Foliage m2/kg as a function
		       ///of (relative) light
  ParametricCurve adf; ///Apical dominance function, apical = f(qin/TreeQinMax)
  ParametricCurve gof; ///Gravelius order effect on segment length 
  ParametricCurve spsd;///Sapwood Down as a function of Gravelius order
  ParametricCurve af;  ///Initial foliage m^2/kgC as function of light
  ParametricCurve naf; ///Needle angle as a function of light
  ParametricCurve wd;  ///Density of growth ring as a function of segment age
  ParametricCurve lr;  ///Segment  length  -  radius relationship  as  a
		       ///function of light: R = f(relative_light)*L
  ParametricCurve ebhf;///Function EBHlambda as a funct. of order
  ParametricCurve bvf; ///Modifier of no. buds as a function of local needle density 

  double hw_start;     ///The age when the heartwood starts to build up,
		       ///defaults to 0
  double hc;           ///The height at the crown limit
  LGMdouble branch_angle; ///Initial angle of the main branch
  double is_EBH;       ///if Extended Borchert &  Honda resource distribution in use
                       ///(0 = false, 1 = true)


};

/////ScotsPineSegment    ///////////////////////

class ScotsPineSegment: public PineSegment<ScotsPineSegment,ScotsPineBud>{
  ///The SetValue  for LGPsf  changes the specific  leaf area to  be a
  ///function instead of being  single tree level parameter. That's why
  ///no  value argument. Also  this is  meant to  be used  with functor
  ///SetScotsPineSegmentSf() only.
  friend LGMdouble SetValue(ScotsPineSegment& ts, LGMAD  name){    //a bit unconventional SetValue
    //The data from P Kaitaniemi suggesta following model for sf (m2/kgC)
    //    sf = 200.0/(5.8307 + 5.3460*relative_height + omega + 27.0*length)
    //Specific  leaf  area depends  oon  the  relative  height of  the
    //segment in a tree, its gravelius order and its length.
    double sf_old = GetValue(ts,name);
    if (name == LGAsf){
      ScotsPineTree& t = dynamic_cast<ScotsPineTree&>(GetTree(ts));
      const ParametricCurve& fsf = GetFunction(t,SPFSF);
      double qin = GetValue(ts,LGAQin);
      double qmax = GetValue(t,TreeQinMax);
      double qrel = qin/qmax;
      double sf = fsf(qrel);
      //cout << "Qin " << qin << " Sf " <<sf << " L " << GetValue(ts,LGAL) <<endl;
      //double t_height = GetValue(t,LGAH);
      //double ts_height = GetValue(ts,LGAH);
      //double relative_height = ts_height/t_height;
      //double omega = GetValue(ts,LGAomega);
      //double length = GetValue(ts,LGAL);
      //Data from P Kaitaniemi for gravelius order, main axis is 1 
      //ParametricCurve forder("1 2.1709 2 1.6257 3 0.4250 4 -0.0724 5 0.000 6 0.00",0);
      //double order = forder(omega);
      //Model(s)  for sf  from  P.  Kaitaniemi data.  The  sf depends  on
      //relative height of the segment, gravelius order and length 
      //double sf = 200.0/(5.8307 + 5.3460*relative_height + order + 27.0*length);
      //double sf = 200.0/(7.3823 + 3.6951*relative_height + order + 27.0*length);
      SetValue(ts,LGAsf,sf);
      return sf_old;
    }
    else{
      cerr << "SetValue(ts,LGAsf) unknown name " << name << endl;
    }
    return sf_old;
  }


  friend LGMdouble GetValue(const ScotsPineSegment& ts, const LGMAD name){
    if (name == LGAWh){
      ///LGAWh: For TreeSegment the sapwood, hertwood and total wood masses are simply rho*LGAV[s,h,wood],
      ///i.e. one density value for all segments. ScotsPineSegment (Lig-Crobas) uses measurement data
      ///for more detailed values. This is a shorthand fix: Hakkila 1969 CIFF 67.6
      ///shows that wood density for branches is > 200 kgC/m3. To achieve this add 10 years to age,
      //if branch
      const ParametricCurve& wd = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(ts)), SPWD);
      double age = GetValue(ts, LGAage);
      double Wh = wd(age)*GetValue(ts,LGAVh);  
      return Wh;
    }   
    else if (name == LGAWood){
      const ParametricCurve& wd = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(ts)), SPWD);
      double age = GetValue(ts, LGAage);
      if (GetValue(ts,LGAomega) > 1.0) {
        age += 10.0;}
      double Wood = wd(age)*GetValue(ts,LGAV);
      return Wood;
    }
    else if (name == LGAWs){
      return GetValue(ts,LGAWood)-GetValue(ts,LGAWh);
    }
    else{
      return GetValue(dynamic_cast<const CfTreeSegment<ScotsPineSegment,ScotsPineBud>&>(ts),name);
    }
  }


  friend LGMdouble GetValue(const ScotsPineSegment& ts, const SPAD  name){
    if (name == SPAAsDown){    
      return ts.AsDown;
    }
    else if (name == SPrue) {
      return ts.rue;
    }

    else{
      cerr << "GetValue(ScotsPineSegment,name ) unknown name " << name <<endl; 
      return  0.0;
    }
  }

  friend LGMdouble SetValue(ScotsPineSegment& ts, const SPAD name, const LGMdouble value)
  {
    LGMdouble old_value = GetValue(ts,name);
    if (name == SPAAsDown)
      ts.AsDown = value;
    else if(name == SPrue) {
      ts.rue = value;
    }
    else{
      cerr << "SetValue(ScotsPineSegment&, SPAD, value) unknown name " << name << endl;
      old_value = 0.0;
    }
    return old_value;
  }
  
public:
  ///\attention LGAQabs initialized to 1.0 due to usage (anything != 0.0 would do)
  ///\sa ScotsPineSegment::photosynthesis
  ScotsPineSegment(const Point& p,const PositionVector& d,
		   const LGMdouble go,const METER l,
		   const METER r,const METER rh,
		   Tree<ScotsPineSegment,ScotsPineBud>* tree)
    :PineSegment<ScotsPineSegment,ScotsPineBud>(p,d,go,l,r,rh,tree),
     EBH_resource(0.0), apical(1.0),rue(1.0)  {SetValue(*this,LGAQabs,1.0);}

  //DiameterGrowth  in   functors  [Try|Do]ScotsPineDiameterGrowth  in
  //DiameterGrowth.h
  TcData& diameterGrowth(TcData& data)
  {
    cerr << "ScotsPineSegment::diameterGrowth no longer in use" <<endl;
    cerr << "Use [Try|Do]ScotsPineDiameterGrowth functors" <<endl;
    cerr << "See DiameterGrowth.h" <<endl;
    return data;
  }

  void respiration();
  void aging();
  void setEBHResource(const double& r) {EBH_resource = r;}
  double getEBHResource() {return EBH_resource;}
  void setQv(const double& qv) {Qv = qv;}
  double getQv() {return Qv;}
  void setApical(const double& ap) {apical = ap;}
  double getApical() {return apical;}
  double view;                   ///Make this private?
  void photosynthesis();
  LGMdouble getQinStand() {return Qin_stand;}
  void setQinStand(LGMdouble qis) {Qin_stand = qis;}
private:
  LGMdouble AsDown;
  LGMdouble Qin_stand;
  LGMdouble EBH_resource;      ///Extended Borchert-Honda model resource
  LGMdouble Qv;      ///Cumulative radiation in basipetal direction in EBH calculation
  LGMdouble apical;  ///Measure of apical dominance on lateral Segments = f(qin/TreeQinMax),
                     ///see functor SetScotsPineSegmentApical
  LGMdouble rue;  ///Radiation use efficiency: photosynthetic production of a CfTreeSegment =
                  ///rue*LGPpr*Qabs, where parameter LGPpr = Photosynthetic efficiency
                  ///(see LGMSymbols.h). rue depends on the radiation conditions of the CfTreeSegment
                  ///at its birth. At full light (at top of the stand) rue = 1, in shaded
                  ///conditions possibly rue > 1.
};   //ScotsPineSegment


class ScotsPineBud:  public PineBud<ScotsPineSegment,ScotsPineBud>{
public:
  ScotsPineBud(const Point& p, const PositionVector& d, 
	       const LGMdouble go, Tree<ScotsPineSegment,ScotsPineBud>* tree)
    :PineBud<ScotsPineSegment,ScotsPineBud>(p,d,go,tree){}
};


class SetScotsPineSegmentLength{
public:
  SetScotsPineSegmentLength(double lamda):l(lamda)  ///space_occupancy is
  {space_occupancy.resetOccupiedTry();}             ///global
  SetScotsPineSegmentLength(const SetScotsPineSegmentLength& sl)
    :l(sl.l)/*,apical(sl.apical),qin(sl.qin)*/{}
  SetScotsPineSegmentLength& operator=(const SetScotsPineSegmentLength& sl){
    l = sl.l;
    // apical = sl.apical;
    // qin = sl.qin;
    return *this;
  }
  TreeCompartment<ScotsPineSegment,ScotsPineBud>* 
  operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0.0){
	double go =  GetValue(*ts,LGAomega); 
	double Lnew = 0.0;

	if(GetValue(dynamic_cast<ScotsPineTree&>(GetTree(*ts)), SPis_EBH) < 1.0) {
	  //The effect of Gravelius order on segment length
	  //	    go = GetValue(*ts,LGAomega); 
	  const ParametricCurve& fgo = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFGO);
	  //Vigour index effect on segment length
	  double vi = GetValue(*ts,LGAvi);
	  //	double q = GetValue(GetTree(*ts),LGPq);
	  //In Tree Physiology for side branches fp is for example as follows:
	  //fp = (1-a)f(vi) = (1-0.2)(0.15+0.85vi) = 0.8(0.15+0.85vi)
	  const ParametricCurve& fvi = GetFunction(GetTree(*ts),LGMVI);
	  //The value of  'apical' is [0,1] for new branches  and set to 1
	  //after that (see below)
	  double my_apical = ts->getApical();
	  // /*double*/ Lnew = apical*fvi(vi)*fgo(go); 
	  Lnew = my_apical*fvi(vi)*fgo(go);
	    
	  /* 	    cout << "a vi fvi go fgo ip " << my_apical << " "<< vi << " "<< fvi(vi) */
	  /* 		 << " "<< go << " "<< fgo(go)<< " "  */
	  /* 		 << GetValue(*ts,LGAQin)/GetValue(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),TreeQinMax) */
	  /* 		 << endl; */
	  //	    double Lnew =  max<double>(0.0,1.0-(go-1.0)*q);
	} else {
	  //This is according to W. Palubicki and K. Horel and S. Longay and
	  //A. Runions and B. Lane and R. Mech and P. Prusinkiewicz. 2009.
	  //Self-organizing tree models for image synthesis ACM Transactions on
	  //Graphics 28 58:1-10. Length growth is according to
	  //Extended Borchert-Honda model, EBH_resource is share [0,1] of total
	  //intercepted radiation that has been distributed according to their function.
	  Lnew = ts->getEBHResource();
	}
      
	//	  cout << "Lnew1 Lnew2 Lnew3 Lnew4 Lnew5  " << Lnew << " ";


	double adhoc_factor = 1.0;
	if(is_adhoc) {
	  if(go > 1.0) {
	    if(global_hcb >= L_H ) {
	      adhoc_factor = 1.0;
	    } else if(global_hcb < 0.0) {
	      adhoc_factor = 1.0;
	    } else {
	      double z = GetPoint(*ts).getZ();
	      double rel_z = (z-global_hcb)/(L_H - global_hcb);
	      if(rel_z > 1.0)
		rel_z = 1.0;
	      if(rel_z < 0.0)
		rel_z = 0.0;
	      adhoc_factor = adhoc(rel_z);
	      if(go < 3.0 && adhoc_factor > 1.5) {
		//  adhoc_factor = 1.0/adhoc_factor;
		adhoc_factor = 1.5;
	      }
	    }
	  }
	}


	//The effect of relative light on segment length
	const ParametricCurve& fip = GetFunction(GetTree(*ts),LGMIP);
	double B = GetValue(GetTree(*ts),TreeQinMax);
	double qin = GetValue(*ts,LGAQin);
	double ip = qin/B;

	
	if(GetValue(dynamic_cast<ScotsPineTree&>(GetTree(*ts)), SPis_EBH) < 1.0) {

	  Lnew = l*fip(ip)*Lnew * adhoc_factor;
	  //	  cout << Lnew << " ";
	} else {

	  Lnew = l * Lnew * adhoc_factor;
	}

	//Random variation in lengths of segments (not stem)
	if(go > 1.0) {
	  extern int ran3_seed;
	  LGMdouble rp = GetValue(GetTree(*ts),LGPlen_random);
	  Lnew *= 1.0 + (rp/0.5)*(ran3(&ran3_seed)-0.5);
	}
	//	  cout << Lnew << " ";

	//Here Space occupancy
	if(space0 || space1 || space2) {
	  if(go > 1.0) {
	    if(Lnew > R_EPSILON && Lnew < 1.5) {
	      Point p0 = GetPoint(*ts);
	      PositionVector dir = GetDirection(*ts);
	      Point end_p = p0 + Lnew*Point(dir);
	      vector<int> e_ind = space_occupancy.getBoxIndexes(end_p);
	      vector<int> s_ind = space_occupancy.getBoxIndexes(p0);

	      if(!((abs(e_ind[0]-s_ind[0])+abs(e_ind[1]-s_ind[1])+abs(e_ind[2]-s_ind[2])) == 0)) {
		bool out_of_space = false;
		VoxelBox& vb = space_occupancy.voxboxes[1][1][1];
		try{
		  vb = space_occupancy.getVoxelBox(end_p);
		}
		catch(Lignum::OutOfVoxelSpaceException& e) {
		  out_of_space = true;
		}

		vector<VoxelBox> neighborhood;
		if(space1) {
		  neighborhood = space_occupancy.
		    getVoxelBoxPositiveNeighborhood(end_p, dir);
		}
		if(space2) {
		  list<vector<int> > neighbors = space_occupancy.
		    getBoxesAroundPoint(end_p, space2_distance, false);
		  list<vector<int> >::iterator I;
		  for(I = neighbors.begin(); I != neighbors.end(); I++) {
		    if(!((abs((*I)[0]-s_ind[0])+abs((*I)[1]-s_ind[1])+abs((*I)[2]-s_ind[2])) == 0)) {
		      neighborhood.push_back(space_occupancy.voxboxes[(*I)[0]][(*I)[1]][(*I)[2]]);
		    }
		  }
		}

		int nb_occupied = 0;
		if(space1 || space2) {
		  int lv = static_cast<int>(neighborhood.size());
		  for(int i = 0; i < lv; i++) {
		    if(neighborhood[i].getOccupied() || neighborhood[i].getOccupiedTry())
		      nb_occupied++;
		  }
		}

		if(!out_of_space) {
		  if(vb.getOccupied() || vb.getOccupiedTry() || (nb_occupied > 0)) {
		    Lnew = 0.0;
		  } else {
		    vb.setOccupiedTry(true);
		  }
		} else {
		  if(nb_occupied > 0) {
		    Lnew = 0.0;
		  }
		}

	      } //if(!((abs(e_ind[0]-s_ind[0])+abs(e_i ...
	    }
	  } //if(go > 1) {
	} //if(space0 || space1 || space2) { 


	Lnew = max(Lnew,0.0);

	const ParametricCurve& flr = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFLR);

	//don't allow segments shorter than LGPmin or very thin, less than 0.2 mm (R = flr * L)
	if (Lnew < GetValue(GetTree(*ts),LGPLmin) || flr(ip)*Lnew < 0.0001) {
	  Lnew = 0.0;
	}
	//	  cout << Lnew << endl;

	SetValue(*ts,LGAL,Lnew);
	//Initial radius
	SetValue(*ts,LGAR,flr(ip)*Lnew);
	//Reset previous Rh!!!!
	SetValue(*ts,LGARh,0.0);
	//Initial heartwood
	SetValue(*ts,LGARh,sqrt((GetValue(GetTree(*ts),LGPxi)*GetValue(*ts,LGAAs))/PI_VALUE));
	//Initial foliage for Scots pine: LGPAf is a function of light
	const ParametricCurve& faf =
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFAF);
	//	  SetValue(*ts,LGAWf,faf(ip)*GetValue(*ts,LGASa));

	SetValue(*ts,LGAWf,faf(ip)*Lnew);   //Needle mass now proportional to segment length

	//Remember the initial foliage!!
	SetValue(*ts,LGAWf0,GetValue(*ts,LGAWf));
	//Remember original sapwood area As0
	SetValue(*ts,LGAAs0,GetValue(*ts,LGAAs));

      } // if (GetValue(*ts,LGAage) == 0.0){

    } //if (ScotsPineSegment* ts ...  i.e. end segment

    return tc;
  }
private:
  double l;//Lamda to iterate segment lengths
};

#endif
