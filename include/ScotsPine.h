#ifndef SCOTSPINE_H
#define SCOTSPINE_H
#include <Pine.h>

//SPFN =  Scots Pine Functions,
//SPFAF = Initial foliage m^2/kgC   
//SPFAD = Scots Pine  Apical Dominance,
//SPFGO = Scots Pine Gravelius Order
//SPFLR = Scots Pine length - radius relationship, R = f(relative_light)*R 
//SPFNA = Scots Pine Needle Angle
//SPFSF = Scots Pine Specific Leaf Area
//SPSD = Scots Pine Sapwood Down as a function of Gravelius order
enum SPFN {SPFAF,SPFAD, SPFGO,SPFLR,SPFNA,SPFSF, SPSD, SPWD};

//Enumeration for SetValue, GetValue in ScotsPine
//Scots Pine Attribute Double SPAD,
//Sapwood down, Height at crown limit, start of heartwood build up 
enum SPAD {SPAAsDown,SPCrownRatio,SPHc,SPHwStart};

// 0  LGAplength  Path length  from the base of the  tree to a segment        
class ScotsPineBud;
class ScotsPineSegment;

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
    default:
      LGMMessage("ScotsPineTree Unknown function");
      throw ParametricCurve();
    }
  }

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
  else{
    cerr << "SetFunction unknown function: " << name << endl;
    throw ParametricCurve();
  }
  return;
}

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
      cerr << "GetValue(ScotsPineTree,name) unknown name " << name <<endl;
    }
    return old_value;
  }
public:
  ScotsPineTree(const Point& p, const PositionVector& d,
		const string& sffun, const string& apicalfun, 
		const string& gofun, const string& sapwdownfun,
		const string& affun, const string& nafun, 
		const string& wdfun, const string& lrfun)
    :Tree<ScotsPineSegment,ScotsPineBud>(p,d),sf(sffun),adf(apicalfun),
     gof(gofun),spsd(sapwdownfun),af(affun),naf(nafun),wd(wdfun),lr(lrfun),
     hw_start(0.0),hc(0.0){}
  LGMdouble getBranchAngle() {return branch_angle;}
  void setBranchAngle(LGMdouble ba) {branch_angle = ba;}

private:
  ParametricCurve sf;//Specific leaf area: Foliage m2/kg as a function
		     //of (relative) light
  ParametricCurve adf;//Apical dominance function, apical = f(qin/TreeQinMax)
  ParametricCurve gof;//Gravelius order effect on segment length 
  ParametricCurve spsd;//Sapwood Down as a function of Gravelius order
  ParametricCurve af;//Initial foliage m^2/kgC as function of light
  ParametricCurve naf;//Needle angle as a function of light
  ParametricCurve wd; //Density of growth ring as a function of segment age
  ParametricCurve lr; //Segment  length  -  radius relationship  as  a
		      //function of light: R = f(relative_light)*L
  double hw_start;    //The age when the heartwood starts to build up,
		      //defaults to 0
  double hc;//The height at the crown limit
  LGMdouble branch_angle;

};

/////ScotsPineSegment    ///////////////////////

class ScotsPineSegment: public PineSegment<ScotsPineSegment,ScotsPineBud>{
  //This SetValue  for LGPsf  changes the specific  leaf area to  be a
  //function instead of being  single tree level parameter. That's why
  //no  value argument. Also  this is  meant to  be used  with functor
  //SetScotsPineSegmentSf() only.
  friend LGMdouble SetValue(ScotsPineSegment& ts,LGMAD  name){
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
      //For TreeSegment the sapwood, hertwood and total wood masses are simply rho*LGAV[s,h,wood],
      //i.e. one density value for all segments. ScotsPineSegment (Lig-Crobas) uses measurement data
      //for more detailed values. This is a shorthand fix: Hakkila 1969 CIFF 67.6
      //shows that wood density for branches is > 200 kgC/m3. To achieve this add 10 years to age,
      //if branch
      const ParametricCurve& wd = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(ts)), SPWD);
      double age = GetValue(ts, LGAage);
      double Wh = wd(age)*GetValue(ts,LGAVh);                                                               
      return Wh;                                                                                        
    }   
    else if (name == LGAWood){                                                                          
      const ParametricCurve& wd = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(ts)), SPWD);         
      double age = GetValue(ts, LGAage);                                                                
      if (GetValue(ts,LGAomega) > 1.0)                                                                  
        age += 10.0;                                                                                    
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
    else
      old_value = 0.0;
    return old_value;
  }


public:
  ScotsPineSegment(const Point& p,const PositionVector& d,
		   const LGMdouble go,const METER l,
		   const METER r,const METER rh,
		   Tree<ScotsPineSegment,ScotsPineBud>* tree)
    :PineSegment<ScotsPineSegment,ScotsPineBud>(p,d,go,l,r,rh,tree)
  {
  }

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
  void photosynthesis();
  LGMdouble getQinStand() {return Qin_stand;}
  void setQinStand(LGMdouble qis) {Qin_stand = qis;}
private:
  LGMdouble AsDown;
  LGMdouble Qin_stand;
};   //ScotsPineSegment






class ScotsPineBud:  public PineBud<ScotsPineSegment,ScotsPineBud>{
public:
  ScotsPineBud(const Point& p, const PositionVector& d, 
	       const LGMdouble go, Tree<ScotsPineSegment,ScotsPineBud>* tree)
    :PineBud<ScotsPineSegment,ScotsPineBud>(p,d,go,tree){}
};


class SetScotsPineSegmentLength{
public:
  SetScotsPineSegmentLength(double lamda):l(lamda),apical(1.0),qin(0.0){}
  SetScotsPineSegmentLength(const SetScotsPineSegmentLength& sl)
    :l(sl.l),apical(sl.apical),qin(sl.qin){}
  SetScotsPineSegmentLength& operator=(const SetScotsPineSegmentLength& sl){
    l = sl.l;
    apical = sl.apical;
    qin = sl.qin;
    return *this;
  }
  TreeCompartment<ScotsPineSegment,ScotsPineBud>* 
  operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0.0){
	const ParametricCurve& fip = GetFunction(GetTree(*ts),LGMIP);
	//The effect of Gravelius order on segment length
	double go = GetValue(*ts,LGAomega); 
	const ParametricCurve& fgo = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFGO);
	//Vigour index effect on segment length
	double vi = GetValue(*ts,LGAvi);
	//double q = GetValue(GetTree(*ts),LGPq);
	//In Tree Physiology for side branches fp is for example as follows:
	//fp = (1-a)f(vi) = (1-0.2)(0.15+0.85vi) = 0.8(0.15+0.85vi)
	const ParametricCurve& fvi = GetFunction(GetTree(*ts),LGMVI);
	//The value of  'apical' is [0,1] for new branches  and set to 1
	//after that (see below)
	double Lnew = apical*fvi(vi)*fgo(go); 
	//double Lnew =  max<double>(0.0,1.0-(go-1.0)*q);
	//The effect of relative light on segment length
	double B = GetValue(GetTree(*ts),TreeQinMax);
	double qin = GetValue(*ts,LGAQin);
	double ip = qin/B;
	
	Lnew = l*fip(ip)*Lnew;

	//Random variation in lengths of segments (not stem)
	if(go > 1.0) {
	  int dummy = 1;
	  LGMdouble rp = GetValue(GetTree(*ts),LGPlen_random);
	  Lnew *= 1.0 + (rp/0.5)*(ran3(&dummy)-0.5);
	}

	Lnew = max(Lnew,0.0);
	//if (go == 1.0)
	//cout << Lnew <<endl;
	//cout << "Lnew " << Lnew << " ip " << ip << " fip(ip) " << fip(ip) << " q " 
	//<< q << " vi " << vi << " fvi " << fvi(vi) << endl;
	//don't allow segments shorter than LGPmin
	if (Lnew < GetValue(GetTree(*ts),LGPLmin))     //////////////////////////////////////
	  Lnew = 0.0;
	SetValue(*ts,LGAL,Lnew);
	//Initial radius
	const ParametricCurve& flr = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFLR);
	SetValue(*ts,LGAR,flr(ip)*Lnew);
	//Reset previous Rh!!!!
	SetValue(*ts,LGARh,0.0);
	//Initial heartwood
	SetValue(*ts,LGARh,sqrt((GetValue(GetTree(*ts),LGPxi)*GetValue(*ts,LGAAs))/PI_VALUE));
	//Initial foliage for Scots pine: LGPAf is a function of light
	const ParametricCurve& faf =
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFAF);
	SetValue(*ts,LGAWf,faf(ip)*GetValue(*ts,LGASa));
	//Remember the initial foliage!!
	SetValue(*ts,LGAWf0,GetValue(*ts,LGAWf));
	//Remember original sapwood area As0
	SetValue(*ts,LGAAs0,GetValue(*ts,LGAAs));
      }
      else{
	//Propagate qin towards newly created segments 
	qin = GetValue(*ts,LGAQin);
      }
    }//end segment
    //Peek  the  segment,  if  it  is new  branching  set  the  apical
    //dominance [0,1] else set it to 1
    else if (Axis<ScotsPineSegment,ScotsPineBud>* axis = 
	     dynamic_cast<Axis<ScotsPineSegment,ScotsPineBud>*>(tc)){
      if(ScotsPineSegment* sps = dynamic_cast<ScotsPineSegment*>(GetFirstTreeCompartment(*axis))) {
	if(sps != NULL && GetValue(*sps, LGAage) == 0.0) {
	  ScotsPineTree& t = dynamic_cast<ScotsPineTree&>(GetTree(*axis));
	  const ParametricCurve& fad = GetFunction(t,SPFAD);
	  //Set the apical dominance as the function of relative light
	  apical = fad(qin/GetValue(GetTree(*axis),TreeQinMax));
	}
      }
      //old branch
      else{
	apical = 1.0;
      }
    }
    return tc;
  }
private:
  double l;//Lamda to iterate segment lengths
  mutable double apical; //Apical dominance, 1 or less, e.g. 0.8
  mutable double qin;//Apical is a function of light
};




#endif
