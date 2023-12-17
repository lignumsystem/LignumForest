#ifndef SCOTSPINE_H
#define SCOTSPINE_H
#include <Pine.h>
#include <VoxelSpace.h>
#include <CrownDensityGlobals.h>
///\file ScotsPine.h
///\brief Implementation of ScotsPineTree, ScotsPineSegment and ScotsPineBud
using namespace PineTree;
using namespace CrownDensity;
using namespace Pine;


namespace LignumForest{
  class ScotsPineTree;
  class ScotsPineSegment;
  class ScotsPineBud;
  ///Scots pine tree
  class ScotsPineTree: public Tree<ScotsPineSegment,ScotsPineBud>{
    ///Return the function asked
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
    ParametricCurve sf;  ///<Specific leaf area: Foliage m2/kg as a function of (relative) light
    ParametricCurve adf; ///<Apical dominance function, apical = f(qin/TreeQinMax)
    ParametricCurve gof; ///<Gravelius order effect on segment length 
    ParametricCurve spsd;///<Sapwood Down as a function of Gravelius order
    ParametricCurve af;  ///<Initial foliage m^2/kgC as function of light
    ParametricCurve naf; ///<Needle angle as a function of light
    ParametricCurve wd;  ///<Density of growth ring as a function of segment age
    ParametricCurve lr;  ///<Segment  length  -  radius relationship  as  a function of light: R = f(relative_light)*L
    ParametricCurve ebhf;///<Function EBHlambda as a funct. of order
    ParametricCurve bvf; ///<Modifier of no. buds as a function of local needle density 
    double hw_start;     ///<The age when the heartwood starts to build up, defaults to 0
    double hc;           ///<The height at the crown limit
    LGMdouble branch_angle; ///<Initial angle of the main branch
    double is_EBH;       ///<If Extended Borchert &  Honda resource distribution in use (0 = false, 1 = true)


  };

  ///ScotsPineSegment    
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
    LGMdouble EBH_resource; ///<Extended Borchert-Honda model resource
    LGMdouble Qv;      ///<Cumulative radiation in basipetal direction in EBH calculation
    LGMdouble apical;  ///<Measure of apical dominance on lateral Segments = f(qin/TreeQinMax),\sa SetScotsPineSegmentApical
    ///Radiation use efficiency: photosynthetic production of a CfTreeSegment =
    ///rue*LGPpr*Qabs, where parameter LGPpr = Photosynthetic efficiency
    ///(see LGMSymbols.h). rue depends on the radiation conditions of the CfTreeSegment
    ///at its birth. At full light (at top of the stand) rue = 1, in shaded
    ///conditions possibly rue > 1.
    LGMdouble rue;
  };   //ScotsPineSegment


  ///ScotsPineBud has no changes compared to PineBud. \sa PineBud
  class ScotsPineBud:  public PineBud<ScotsPineSegment,ScotsPineBud>{
  public:
    ScotsPineBud(const Point& p, const PositionVector& d, 
		 const LGMdouble go, Tree<ScotsPineSegment,ScotsPineBud>* tree)
      :PineBud<ScotsPineSegment,ScotsPineBud>(p,d,go,tree){}
  };


  ///\brief Adjust lengths \f$ L \f$ for newly created segments.
  ///
  ///\test Segment length models to allocate resources.
  ///Choose and adapt appropriate models with command line options.<br>
  ///Currently the functor captures the following experiments.
  ///-# Basic model where \f$ L = \lambda \times f_{ip} \times f_{go} \times f_{vi} \f$
  ///-# So called EBH model to distribute resources. Mutually exclusive with Basic model.
  ///-# Optional Ad hoc submodel (ParametricCurve) based on relative height to have variability
  ///   in segment lengths (even for segments with similar conditions)
  ///-# Optional space occupancy submodel to set segment length 0  if segment in already occupied space
  ///-# Optional hard coded random effect in branches (both Basic and EBH model)
  ///-# Growth mode and architecture mode change in Basic model
  ///    -# Architecture mode change in L-system  tries to generate flat branches arranged in a plane
  ///    -# Growth mode change applies new  \f$f_{ip}\f$, \f$f_{go}\f$ functions and tree parameters (*Tree.txt*)
  ///    -# Growth mode and architecture mode are independent from each other.
  ///\sa CrownDensity::Usage 
  class SetScotsPineSegmentLength{
  public:
    ///\deprecated fgo and fip should be set explicitely
    ///\param lamda Lambda to iterate segment length
    SetScotsPineSegmentLength(double lamda):l(lamda)  
    {space_occupancy.resetOccupiedTry();}
    ///Constructor setting Gravelius order (fgo) and relative light (fip) functions
    ///to model growth mode change.
    ///\param lamda Lambda to iterate segment length
    ///\param fgo1 Gravelius order function on segment length
    ///\param fip1 Relative light function on sgment length
    ///\param fgo2 Gravelius order function after growth mode change
    ///\param fip2 Relative light function  after growth mode change
    ///\note fgo1, fip1, fgo2 and fip2 can be set in the main loop   
    SetScotsPineSegmentLength(double lamda,const ParametricCurve& fgo1, const ParametricCurve& fip1,
			      const ParametricCurve& fgo2, const ParametricCurve& fip2)
      :l(lamda),fgo(fgo1),fip(fip1),fgo_mode(fgo2),fip_mode(fip2)  
    {space_occupancy.resetOccupiedTry();}
    SetScotsPineSegmentLength(const SetScotsPineSegmentLength& sl)
      :l(sl.l),fgo(sl.fgo),fip(sl.fip),fgo_mode(sl.fgo_mode),fip_mode(sl.fip_mode)
       /*,apical(sl.apical),qin(sl.qin)*/{}
    SetScotsPineSegmentLength& operator=(const SetScotsPineSegmentLength& sl){
      l = sl.l;
      fgo=sl.fgo;
      fip = sl.fip;
      fgo_mode = sl.fgo_mode;
      fip_mode = sl.fip_mode;
      // apical = sl.apical;
      // qin = sl.qin;
      return *this;
    }
    ///\brief Set segment new length.
    ///
    ///Experiment with various methods and mechanisms to set segment length. Set appropriate command line options.
    ///\test Basic model  \f$L = \lambda \times R \times A \times f_{ip} \times f_{go} \times f_{vi}\f$
    ///where \f$R\f$ is the random effect and \f$A\f$ is apical dominance. \sa ScotsPineSegment::getApical LGPlen_random for a tree
    ///\test EBH model
    ///\test Space occupancy model (segment length on/off)
    ///\test Growth mode change based on Basic model
    ///\param tc Tree compartment
    ///\note `CrownDensity::is_mode_change` and `CrownDensity::mode_change_year` are global variables
    ///\sa CrownDensity::is_architecture_change CrownDensity::architecture_change_year
    ///\sa CrownDensity::Usage
    TreeCompartment<ScotsPineSegment,ScotsPineBud>* 
    operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
    { ///\section setsegmentlength Steps in setting segment length 
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	if (GetValue(*ts,LGAage) == 0.0){
	  double go =  GetValue(*ts,LGAomega); 
	  double Lnew = 0.0;
	  ///\par Set growth functions
	  ///Functions before and after growth mode change, fip(ip) and fgo(go).
	  ///\snippet{lineno} ScotsPine.h SetFuncs
	  ///\internal
	  // [SetFuncs]
	  //First set default f(go) and f(vi), no growth mode chane yet
	  ParametricCurve fgo_tmp(fgo);
	  ParametricCurve fip_tmp(fip);
	  //If growth mode change has happened set appropriate f(go) and f(vi)
	  if (is_mode_change && GetValue(GetTree(*ts),LGAage) >= mode_change_year){
	    fgo_tmp = fgo_mode;
	    fip_tmp = fip_mode;
	  }
	  // [SetFuncs]
	  ///\endinternal
	  //No EBH -> Basic model
	  if(GetValue(dynamic_cast<ScotsPineTree&>(GetTree(*ts)), SPis_EBH) < 1.0) {
	    //Vigour index effect on segment length
	    double vi = GetValue(*ts,LGAvi);
	    //	double q = GetValue(GetTree(*ts),LGPq);
	    //In Tree Physiology for side branches fp is for example as follows:
	    //fp = (1-a)f(vi) = (1-0.2)(0.15+0.85vi) = 0.8(0.15+0.85vi)
	    const ParametricCurve& fvi = GetFunction(GetTree(*ts),LGMVI);
	    //The value of  'apical' is [0,1] for new branches  and set to 1
	    //after that (see below)
	    ///\par Basic segment elongation model 
	    ///\snippet{lineno} ScotsPine.h LBasic1
	    ///\internal
	    // [LBasic1]
	    double my_apical = ts->getApical();
	    //Intermediate result for Basic model (no EBH)
	    //Apical, vigour index and Gravelius order
	    Lnew = my_apical*fvi(vi)*fgo_tmp(go);
	    // [LBasic1]
	    ///\endinternal
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
	    //NOTE: resers/overrides previous LNew
	    ///\par EBH segment elongation model
	    ///\snippet{lineno} ScotsPine.h EBH1
	    ///\internal
	    // [EBH1]
	    // EBH model in use (no Basic model)
	    Lnew = ts->getEBHResource();
	    // [EBH1]
	    ///\endinternal
	  }
	  //	  cout << "Lnew1 Lnew2 Lnew3 Lnew4 Lnew5  " << Lnew << " ";
	  ///\par Ad hoc addition segment elongation model 
	  ///\snippet{lineno} ScotsPine.h ADHOC
	  ///\internal
	  // [ADHOC]
	  //Calculate adhoc factor for segment length i branches 
	  double adhoc_factor = 1.0;
	  if(is_adhoc) {
	    if(go > 1.0) {
	      //Check bounding values for adhoc factor
	      if(global_hcb >= L_H ) {
		adhoc_factor = 1.0;
	      } else if(global_hcb < 0.0) {
		adhoc_factor = 1.0;
	      } else {
		//Find relative heigth value
		double z = GetPoint(*ts).getZ();
		double rel_z = (z-global_hcb)/(L_H - global_hcb);
		if(rel_z > 1.0)
		  rel_z = 1.0;
		if(rel_z < 0.0)
		  rel_z = 0.0;
		//adhoc is function of the relative heigth value
		adhoc_factor = adhoc(rel_z);
		//Final adjustment for adhoc value
		if(go < 3.0 && adhoc_factor > 1.5) {
		  //  adhoc_factor = 1.0/adhoc_factor;
		  adhoc_factor = 1.5;
		}
	      }
	    }
	  }
	  // [ADHOC]
	  ///\endinternal
	  //Relative light
	  double B = GetValue(GetTree(*ts),TreeQinMax);
	  double qin = GetValue(*ts,LGAQin);
	  double ip = qin/B;
	  //If basic model (no EBH) add f(ip) and  *ad_hoc* to Lnew
	  if(GetValue(dynamic_cast<ScotsPineTree&>(GetTree(*ts)), SPis_EBH) < 1.0) {
	    ///\par Relative light addition to Basic segment elongation
	    ///Default value for `adhoc_factor` is 1
	    ///\snippet{lineno} ScotsPine.h LBasic2
	    ///\internal
	    // [LBasic2]
	    //Result of the Basic model if no ad_hoc
	    //The l denotes lambda, fip(ip) denotes relative light
	    Lnew = l*fip_tmp(ip)*adhoc_factor*Lnew;
	    // [LBasic2]
	    ///\endinternal
	    //	  cout << Lnew << " ";
	  } else {
	    ///\par Ad hoc addition to EBH  segment elongation
	    ///\note No relative light condition with fip(ip)
	    ///\snippet{lineno} ScotsPine.h EBH2
	    ///\internal
	    // [EBH2]
	    //Result for EBH, no f(ip) only ad_hoc
	    //The l denotes lambda
	    Lnew = l*adhoc_factor*Lnew;
	    // [EBH2]
	    ///\endinternal
	  }

	  //Random variation in lengths of segments (not stem)
	  if(go > 1.0) {
	    LGMdouble rp = GetValue(GetTree(*ts),LGPlen_random);
	    ///\par Optional random segment length variation
	    ///Applied to both Basic and EBH models
	    ///\snippet{lineno} ScotsPine.h RANDOM
	    ///\internal
	    // [RANDOM]
	    // For branches if random variation in use
	    Lnew *= 1.0 + (rp/0.5)*(ran3(&CrownDensity::ran3_seed)-0.5);
	    // [RANDOM]
	    ///\endinternal
	  }

	  //Here Space occupancy, segment length on/off
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
		  ///\par Optional space occupancy check
		  ///Set segment length to 0 if growth space (voxel) occupied
		  ///\snippet{lineno} ScotsPine.h SPACE
		  ///\internal
		  // [SPACE]
		  //If growth space checking is on
		  //the segment length L is set to 0
		  //if the growth space is occupied
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
		  // [SPACE]
		  ///\endinternal
		} //if(!((abs(e_ind[0]-s_ind[0])+abs(e_i ...
	      }
	    } //if(go > 1) {
	  } //if(space0 || space1 || space2) { 
	  ///\par Segment length consitency checks
	  ///\snippet{lineno} ScotsPine.h LCHECKS
	  ///\internal
	  // [LCHECKS]
	  // Final checks
	  //1. No negative lengths
	  Lnew = max(Lnew,0.0);
	  //NOTE: Segment length/radius function as function of relative light
	  const ParametricCurve& flr = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFLR);
	  //2. Don't allow segments shorter than LGPmin or very thin segments, less than 0.2 mm (R = flr * L)
	  if (Lnew < GetValue(GetTree(*ts),LGPLmin) || flr(ip)*Lnew < 0.0001) {
	    Lnew = 0.0;
	  }
	  // [LCHECKS]
	  ///\endinternal
	  //	  cout << Lnew << endl;
	  ///\par Set segment dimensions
	  ///\snippet{lineno} ScotsPine.h SDIMENSIONS
	  ///\internal
	  // [SDIMENSIONS]
	  //Set segment dimensions based on segment length
	  SetValue(*ts,LGAL,Lnew);
	  //Initial radius, NOTE: R = flr(ip)*Lnew, i,.e. length/radius function NOT static constant  
	  SetValue(*ts,LGAR,flr(ip)*Lnew);
	  //Reset previous Rh!!!!
	  SetValue(*ts,LGARh,0.0);
	  //Initial heartwood
	  SetValue(*ts,LGARh,sqrt((GetValue(GetTree(*ts),LGPxi)*GetValue(*ts,LGAAs))/PI_VALUE));
	  //Initial foliage for Scots pine: LGPAf is a function of light
	  const ParametricCurve& faf = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFAF);
	  //	  SetValue(*ts,LGAWf,faf(ip)*GetValue(*ts,LGASa));
	  //Needle mass now proportional to segment length
	  SetValue(*ts,LGAWf,faf(ip)*Lnew);   
	  //Remember the initial foliage!!
	  SetValue(*ts,LGAWf0,GetValue(*ts,LGAWf));
	  //Remember original sapwood area As0
	  SetValue(*ts,LGAAs0,GetValue(*ts,LGAAs));
	  // [SDIMENSIONS]
	  ///\endinternal
	} // if (GetValue(*ts,LGAage) == 0.0){

      } //if (ScotsPineSegment* ts ...  i.e. end segment

      return tc;
    }
  private:
    double l;///<Lamda to iterate segment lengths
    ParametricCurve fgo; ///< Gravelius order function
    ParametricCurve fip; ///< fip (rel. light) function
    ParametricCurve fgo_mode;///< Gravelius order function after growth mode change
    ParametricCurve fip_mode;///< fip (rel. light) function after growth mode change
  };


  class FindMaxRueQin{
  public:
    LGMdouble&
    operator()(LGMdouble& max_rue, TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const 
    {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	LGMdouble qin = GetValue(*ts,LGAQin);
	LGMdouble rue = GetValue(*ts,SPrue);
	if(rue * qin > max_rue) {
	  max_rue = qin * rue;
	}
      }
      return max_rue;
    }
  };

}//End namespace LignumForest
#endif
