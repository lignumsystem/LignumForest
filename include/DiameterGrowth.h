#ifndef DIAMETERGROWTH_H
#define DIAMETERGROWTH_H
#include <Lignum.h>
#include <ScotsPine.h>

// PartialSapwoodAreaDown makes it possible  to pass sapwood down based
// on the  gravelius order of  the segments. The percentage  of sapwood
// down is defined in the ParametricCurve given in the constructor.

// Usage: AccumulateDown(tree,DiameterGrowthData,PartialSapwoodAreaDown(ParametricCurve),DiameterGrowth())

// PartialSapwoodAreaDown is the user defined "add and assign" operator
// that is  called by  AccumulateDown at each  branching point.   d1 is
// coming downwards from the segment s1 above branching point b1 in the
// same axis, i.e.   s1 and the b1 have the  same gravelius order.  The
// 'd2' is coming  from the last segment  s2 in the axis in  b1, i.e s2
// has one  order higher gravelius order  than b1.  d2 is  summed to d1
// and d1 is passed downwards to  the segment s2 below in the same axis
// (if s1 is the last one in an axis, d1 is passed to a branching point
// the axis  belongs to  as d2).  This  summation is repeated  for each
// axis in the branching point.  The functor DiameterGrowth() should at
// each Bud initialize the gravelius order and at each segment s3 below
// b1 implement the diameter growth below,  i.e. s1, b1 and s3 have the
// same gravelius order
class PartialSapwoodAreaDown{
public:
  PartialSapwoodAreaDown(const PartialSapwoodAreaDown& swdown):fsapwdown(swdown.fsapwdown){}
  PartialSapwoodAreaDown(const ParametricCurve& f):fsapwdown(f){}
  DiameterGrowthData& operator()(DiameterGrowthData& d1, DiameterGrowthData& d2)const
  {
    double o1 = GetValue(d1,LGAomega);
    double o2 = GetValue(d2,LGAomega);
    if (o1 < o2){
      double percent = fsapwdown(o2);
      double As = GetValue(d2,LGAAs);
      double Asdown = percent*As;
      SetValue(d1,LGAAs,GetValue(d1,LGAAs)+Asdown);
      SetValue(d1,DGWs,GetValue(d1,DGWs)+GetValue(d2,DGWs));
      SetValue(d1,DGWf,GetValue(d1,DGWf)+GetValue(d2,DGWf));
      SetValue(d1,DGWfnew,GetValue(d1,DGWfnew)+GetValue(d2,DGWfnew));
    }
    else{
      //cout << "PartialSapwoodAreaDown " << o1 << " >= " << o2 <<endl;
      //cout << "Passing all spawood down " << GetValue(d2,LGAAs) << " " << GetValue(d2,DGWf) << " " 
      //   << GetValue(d2,DGWfnew) << endl;
      SetValue(d1,LGAAs,GetValue(d1,LGAAs)+GetValue(d2,LGAAs));
      SetValue(d1,DGWs,GetValue(d1,DGWs)+GetValue(d2,DGWs));
      SetValue(d1,DGWf,GetValue(d1,DGWf)+GetValue(d2,DGWf));
      SetValue(d1,DGWfnew,GetValue(d1,DGWfnew)+GetValue(d2,DGWfnew));
    }
    //cout << "O1 " << GetValue(d1,LGAomega) <<endl;
    return d1;
  }
private:
  const ParametricCurve& fsapwdown;
};

//This is must be the same as DoScotsPineDiameterGrowth method, but we
//can't change the segment's dimensions.

class TryScotsPineDiameterGrowth{
public:
  DiameterGrowthData& operator()(DiameterGrowthData& data, TreeCompartment<ScotsPineSegment,
				 ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0){//New segment
	const ParametricCurve& wd = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPWD);
	double rho = wd(0.0);
	double newWs = rho*GetValue(*ts, LGAAs)*GetValue(*ts,LGAL);
        //Sapwood down as a function of gravelius order (LGAomega)                                                
        //Main axis is 1, first branches 2 etc.                                                                   
        const ParametricCurve& sapwdown = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPSD);           
        //Berninger and Nikinmaa Can J suggests that about 0.75 of the
        //sapwood area goes  down when LGAomega = 2,  Nikinmaa pers com
        //suggests that a bit more, 0.89, goes down in subbranches All
        //goes down in the main axis
        SetValue(data,LGAAs,
		 GetValue(*ts,LGAAs)*sapwdown(GetValue(*ts,LGAomega)));
        //Collect masses to account for in P-M=G
	//New foliage is a cost
        SetValue(data,DGWfnew,GetValue(*ts,LGAWf));
	//Existing (total)  foliage is for bookkeeping
        SetValue(data,DGWf,GetValue(*ts,LGAWf));
        SetValue(data,DGWs,newWs);
      }
      else{//old segment
        const ParametricCurve& fm = GetFunction(GetTree(*ts),LGMFM);
        //Sapwood requirement of  remaining foliage, assume fm returns
        //proportion of initial foliage present, declining function of
        //age from 1 to 0.
        LGMdouble Asr = fm(GetValue(*ts,LGAage))*GetValue(*ts,LGAAs0);
        //sapwood area from above
        LGMdouble Asu = GetValue(data,LGAAs); 
        //own heartwood, assume aging has done
        LGMdouble Ahown  = GetValue(*ts,LGAAh);
        //requirement for new radius:  sapwood above + own heartwood +
        //own foliage
        LGMdouble Rnew = sqrt((Asu + Ahown + Asr)/PI_VALUE);
        //compare Rnew to R, choose max
        Rnew = max(Rnew, GetValue(*ts,LGAR));
        //New sapwood requirement, thickness growth
        double Asnew = PI_VALUE*pow(Rnew,2.0) -  GetValue(*ts,LGAA);

        //Mass of  the new  sapwood This is  a shorthand  fix: Hakkila
	//1969 CIFF 67.6 shows that wood density for branches is > 200
	//kgC/m3. To achieve this add 10 years to age, if branch
	const ParametricCurve& wd = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPWD);
	double age = GetValue(*tc, LGAage);
	if (GetValue(*ts,LGAomega) > 1.0)
	  age += 10.0;	
	double rho = wd(age);
	double Wsnew = rho*Asnew*GetValue(*ts,LGAL);

        //Sapwood down as a function of gravelius order (LGAomega)
        //Main axis is 1, first branches 2 etc.
        const ParametricCurve& sapwdown = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)),SPSD);
        //Berninger and Nikinmaa Can J suggests that about 0.75 of the
        //sapwood area goes down when LGAomeaga = 2, Nikinmaa pers com
        //suggests that a bit more, 0.89, goes down in subbranches All
        //goes down in the main axis
        SetValue(data,LGAAs,
		 (Asnew+GetValue(*ts,LGAAs))*sapwdown(GetValue(*ts,LGAomega)));
        //Mass of sapwood used in diamater growth
        SetValue(data,DGWs,GetValue(data,DGWs)+Wsnew);
        //Total foliage
        SetValue(data,DGWf,GetValue(data,DGWf)+GetValue(*ts,LGAWf));
      }
    }
    return data;
  }
};

class DoScotsPineDiameterGrowth{
public:
  DiameterGrowthData& operator()(DiameterGrowthData& data,TreeCompartment<ScotsPineSegment,
				 ScotsPineBud>* tc)const 
  {
    //Dimensions of the new segments (age == 0) are iteratively set.
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0.0){
	//Mass of the new sapwood
	double Wsnew = GetValue(*ts,LGAWs);
	//Sapwood down as a function of gravelius order (LGAomega)
	//Main axis is 1, first branches 2 etc.
	const ParametricCurve& sapwdown = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)), SPSD);
	//Berninger and Nikinmaa Can J suggests that about 0.75 of the
	//sapwood area goes down when  LGAomega = 2, Nikinmaa pers com
	//suggests that a bit more, 0.89, goes down in subbranches All
	//goes down in the main axis
	SetValue(data,LGAAs,GetValue(*ts,LGAAs)*sapwdown(GetValue(*ts,LGAomega)));
	//Collect masses for bookkeeping
	SetValue(data,DGWfnew,GetValue(*ts,LGAWf));
	SetValue(data,DGWf,GetValue(*ts,LGAWf));
	SetValue(data,DGWs,Wsnew);
      }
      //Set the dimensions of the existing segments
      else{// (GetValue(*ts,LGAage) > 0.0)
	const ParametricCurve& fm = GetFunction(GetTree(*ts),LGMFM);
	//Sapwood requirement of  remaining foliage, assume fm returns
	//proportion initial foliage  present, declining function from
	//1 to 0.
	LGMdouble Asr = fm(GetValue(*ts,LGAage))*GetValue(*ts,LGAAs0);
	LGMdouble Asu = GetValue(data,LGAAs); //sapwood area from above
	LGMdouble Ahown  = GetValue(*ts,LGAAh);//own heartwood
	//possible new radius
	LGMdouble Rnew = sqrt((Asu + Ahown + Asr)/PI_VALUE);
	//compare Rnew to R, choose max
	Rnew = max(Rnew, GetValue(*ts,LGAR));
	//New sapwood requirement, thickness growth
	double Asnew = PI_VALUE*pow(Rnew,2.0) -  GetValue(*ts,LGAA);
	//Set segment radius
	SetValue(*ts,LGAR,Rnew);
	//Mass of  the new  sapwood This is  a shorthand  fix: Hakkila
	//1969 CIFF 67.6 shows that wood density for branches is > 200
	//kgC/m3. To achieve this add 10 years to age, if branch
	const ParametricCurve& wd = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPWD);
	double age = GetValue(*tc, LGAage);
	if (GetValue(*ts,LGAomega) > 1.0)
	  age += 10.0;	
	double rho = wd(age);
	double Wsnew = rho*Asnew*GetValue(*ts,LGAL);
	//Sapwood down as a function of gravelius order (LGAomega)
	//Main axis is 1, first branches 2 etc.
	const ParametricCurve& sapwdown = 
	  GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)), SPSD);
	//Berninger and Nikinmaa Can J suggests that about 0.75 of the
	//sapwood area goes down when LGAomeaga = 2, Nikinmaa pers com
	//suggests that a bit more, 0.89, goes down in subbranches All
	//goes down in the main axis
	SetValue(data,LGAAs,GetValue(*ts,LGAAs)*sapwdown(GetValue(*ts,LGAomega)));
	//Mass of sapwood used in diamater growth
	SetValue(data,DGWs,GetValue(data,DGWs)+Wsnew);
	//Total foliage
	SetValue(data,DGWf,GetValue(data,DGWf)+GetValue(*ts,LGAWf));
      }
    }
    return data;
  }
};

//This    ScotsPineDiameterGrowth   functor    can   be    used   with
//LGMGrowthAllocator    implemented    in   LGMGrowthAllocator.h    in
//stl-lignum. The  idea is to have  the two modes,  allocation and the
//growth  itself, together  to make  the allocation  of photosynthates
//less error prone
class ScotsPineDiameterGrowth{
public:
  ScotsPineDiameterGrowth(LGMALLOCATORMODE m):mode(m){}
  DiameterGrowthData& operator()(DiameterGrowthData& data,TreeCompartment<ScotsPineSegment,
				 ScotsPineBud>* tc)const
  {
    if (mode == LGMALLOCATE){
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	if (GetValue(*ts,LGAage) == 0){//New segment
	  const ParametricCurve& wd = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPWD);
	  double rho = wd(0.0);
	  double newWs = rho*GetValue(*ts, LGAAs)*GetValue(*ts,LGAL);
	  //Sapwood down as a function of gravelius order (LGAomega)                                     
	  //Main axis is 1, first branches 2 etc.                                                        
	  const ParametricCurve& sapwdown = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPSD);           
	  //Berninger and Nikinmaa Can J suggests that about 0.75 of the
	  //sapwood area goes  down when LGAomega = 2,  Nikinmaa pers com
	  //suggests that a bit more, 0.89, goes down in subbranches All
	  //goes down in the main axis
	  SetValue(data,LGAAs,
		   GetValue(*ts,LGAAs)*sapwdown(GetValue(*ts,LGAomega)));
	  //Collect masses to account for in P-M=G
	  //New foliage is a cost
	  SetValue(data,LGAiWf,GetValue(*ts,LGAWf));
	  //Existing (total)  foliage is for bookkeeping
	  SetValue(data,DGWf,GetValue(*ts,LGAWf));
	  SetValue(data,LGAiWs,newWs);
	}
	else{//old segment
	  const ParametricCurve& fm = GetFunction(GetTree(*ts),LGMFM);
	  //Sapwood requirement of  remaining foliage, assume fm returns
	  //proportion of initial foliage present, declining function of
	  //age from 1 to 0.
	  LGMdouble Asr = fm(GetValue(*ts,LGAage))*GetValue(*ts,LGAAs0);
	  //sapwood area from above
	  LGMdouble Asu = GetValue(data,LGAAs); 
	  //own heartwood, assume aging has done
	  LGMdouble Ahown  = GetValue(*ts,LGAAh);
	  //requirement for new radius:  sapwood above + own heartwood +
	  //own foliage
	  LGMdouble Rnew = sqrt((Asu + Ahown + Asr)/PI_VALUE);
	  //compare Rnew to R, choose max
	  Rnew = max(Rnew, GetValue(*ts,LGAR));
	  //New sapwood requirement, thickness growth
	  double Asnew = PI_VALUE*pow(Rnew,2.0) -  GetValue(*ts,LGAA);
	  //Mass of  the new  sapwood This is  a shorthand  fix: Hakkila
	  //1969 CIFF 67.6 shows that wood density for branches is > 200
	  //kgC/m3. To achieve this add 10 years to age, if branch
	  const ParametricCurve& wd = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPWD);
	  double age = GetValue(*tc, LGAage);
	  if (GetValue(*ts,LGAomega) > 1.0)
	    age += 10.0;	
	  double rho = wd(age);
	  double Wsnew = rho*Asnew*GetValue(*ts,LGAL);
	  //Sapwood down as a function of gravelius order (LGAomega)
	  //Main axis is 1, first branches 2 etc.
	  const ParametricCurve& sapwdown = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)),SPSD);
	  //Berninger and Nikinmaa Can J suggests that about 0.75 of the
	  //sapwood area goes down when LGAomeaga = 2, Nikinmaa pers com
	  //suggests that a bit more, 0.89, goes down in subbranches All
	  //goes down in the main axis
	  SetValue(data,LGAAs,
		   (Asnew+GetValue(*ts,LGAAs))*sapwdown(GetValue(*ts,LGAomega)));
	  //Mass of sapwood used in diamater growth
	  SetValue(data,LGAiWs,GetValue(data,LGAiWs)+Wsnew);
	  //Total foliage
	  SetValue(data,DGWf,GetValue(data,DGWf)+GetValue(*ts,LGAWf));
	}
      }
      return data;
    }
    else {//(mode == LGMGROWTH)
      //Dimensions of the new segments (age == 0) are iteratively set.
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	if (GetValue(*ts,LGAage) == 0.0){
	  //Mass of the new sapwood
	  double Wsnew = GetValue(*ts,LGAWs);
	  //Sapwood down as a function of gravelius order (LGAomega)
	  //Main axis is 1, first branches 2 etc.
	  const ParametricCurve& sapwdown = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)), SPSD);
	  //Berninger and Nikinmaa Can J suggests that about 0.75 of the
	  //sapwood area goes down when  LGAomega = 2, Nikinmaa pers com
	  //suggests that a bit more, 0.89, goes down in subbranches All
	  //goes down in the main axis
	  SetValue(data,LGAAs,GetValue(*ts,LGAAs)*sapwdown(GetValue(*ts,LGAomega)));
	  //Collect masses for bookkeeping
	  SetValue(data,LGAiWf,GetValue(*ts,LGAWf));
	  SetValue(data,DGWf,GetValue(*ts,LGAWf));
	  SetValue(data,LGAiWs,Wsnew);
	}
	//Set the dimensions of the existing segments
	else{// (GetValue(*ts,LGAage) > 0.0)
	  const ParametricCurve& fm = GetFunction(GetTree(*ts),LGMFM);
	  //Sapwood requirement of  remaining foliage, assume fm returns
	  //proportion initial foliage  present, declining function from
	  //1 to 0.
	  LGMdouble Asr = fm(GetValue(*ts,LGAage))*GetValue(*ts,LGAAs0);
	  LGMdouble Asu = GetValue(data,LGAAs); //sapwood area from above
	  LGMdouble Ahown  = GetValue(*ts,LGAAh);//own heartwood
	  //possible new radius
	  LGMdouble Rnew = sqrt((Asu + Ahown + Asr)/PI_VALUE);
	  //compare Rnew to R, choose max
	  Rnew = max(Rnew, GetValue(*ts,LGAR));
	  //New sapwood requirement, thickness growth
	  double Asnew = PI_VALUE*pow(Rnew,2.0) -  GetValue(*ts,LGAA);
	  //New wood radius
	  SetValue(*ts,LGAR,Rnew);
	  //Mass of  the new  sapwood This is  a shorthand  fix: Hakkila
	  //1969 CIFF 67.6 shows that wood density for branches is > 200
	  //kgC/m3. To achieve this add 10 years to age, if branch
	  const ParametricCurve& wd = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPWD);
	  double age = GetValue(*tc, LGAage);
	  if (GetValue(*ts,LGAomega) > 1.0)
	    age += 10.0;	
	  double rho = wd(age);
	  double Wsnew = rho*Asnew*GetValue(*ts,LGAL);
	  //Sapwood down as a function of gravelius order (LGAomega)
	  //Main axis is 1, first branches 2 etc.
	  const ParametricCurve& sapwdown = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)), SPSD);
	  //Berninger and Nikinmaa Can J suggests that about 0.75 of the
	  //sapwood area goes down when LGAomeaga = 2, Nikinmaa pers com
	  //suggests that a bit more, 0.89, goes down in subbranches All
	  //goes down in the main axis
	  SetValue(data,LGAAs,GetValue(*ts,LGAAs)*sapwdown(GetValue(*ts,LGAomega)));
	  //Mass of sapwood used in diamater growth
	  SetValue(data,LGAiWs,GetValue(data,LGAiWs)+Wsnew);
	  //Total foliage
	  SetValue(data,DGWf,GetValue(data,DGWf)+GetValue(*ts,LGAWf));
	}
      }
      return data;
    }
  }
private:
  LGMALLOCATORMODE mode;
};


// This   ScotsPineDiameterGrowth2    functor   can   be    used   with
// PartialSapwoodAreaDown   implemented   in    this   file   and   the
// LGMGrowthAllocator2 implemented in  stl-lignum.  This functor differs
// from ScotsPineDiameterGrowth  in a sense that it  passes all sapwood
// area  down  between  segments. PartialSapwoodAreaDown  cuts  sapwood
// requirement in  branching points from the segments  that have higher
// gravelius order. The result is  that all sapwood area is passed down
// between segments in the same axis.

// Note that there are  two LGMGrowthAllocator functors. The first one
// can  be  used   if  user  relies  on  '+='   operator  defined  for
// AccumulateDown. The  second one accepts user  defined functor (like
// PartialSapwoodAreaDown).  For the details  see LGMGrowthAllocator.h
// in stl-lignum.
class ScotsPineDiameterGrowth2{
public:
  ScotsPineDiameterGrowth2(LGMALLOCATORMODE m):mode(m){}
  DiameterGrowthData& operator()(DiameterGrowthData& data,TreeCompartment<ScotsPineSegment,
				 ScotsPineBud>* tc)const
  {
    //Set gravelius order to be  used at branching points to determine
    //how much sapwood is required from segments above
    if (ScotsPineBud* b = dynamic_cast<ScotsPineBud*>(tc)){
      SetValue(data,LGAomega,GetValue(*b,LGAomega));
      return data;
    }
    if (mode == LGMALLOCATE){
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	if (GetValue(*ts,LGAage) == 0){//New segment
	  const ParametricCurve& wd = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPWD);
	  double rho = wd(0.0);
	  double newWs = rho*GetValue(*ts, LGAAs)*GetValue(*ts,LGAL);
	  SetValue(data,LGAAs,GetValue(*ts,LGAAs));
	  //Collect masses to account for in P-M=G
	  //New foliage is a cost
	  SetValue(data,LGAiWf,GetValue(*ts,LGAWf));
	  //Existing (total)  foliage is for bookkeeping
	  SetValue(data,DGWf,GetValue(*ts,LGAWf));
	  SetValue(data,LGAiWs,newWs);
	}
	else{//old segment
	  const ParametricCurve& fm = GetFunction(GetTree(*ts),LGMFM);
	  //Sapwood requirement of  remaining foliage, assume fm returns
	  //proportion of initial foliage present, declining function of
	  //age from 1 to 0.
	  LGMdouble Asr = fm(GetValue(*ts,LGAage))*GetValue(*ts,LGAAs0);
	  //sapwood area from above
	  LGMdouble Asu = GetValue(data,LGAAs); 
	  //own heartwood, assume aging has done
	  LGMdouble Ahown  = GetValue(*ts,LGAAh);
	  //requirement for new radius:  sapwood above + own heartwood +
	  //own foliage
	  LGMdouble Rnew = sqrt((Asu + Ahown + Asr)/PI_VALUE);
	  //compare Rnew to R, choose max
	  Rnew = max(Rnew, GetValue(*ts,LGAR));
	  //New sapwood requirement, thickness growth
	  double Asnew = PI_VALUE*pow(Rnew,2.0) -  GetValue(*ts,LGAA);
	  //Mass of  the new  sapwood This is  a shorthand  fix: Hakkila
	  //1969 CIFF 67.6 shows that wood density for branches is > 200
	  //kgC/m3. To achieve this add 10 years to age, if branch
	  const ParametricCurve& wd = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPWD);
	  double age = GetValue(*tc, LGAage);
	  if (GetValue(*ts,LGAomega) > 1.0)
	    age += 10.0;	
	  double rho = wd(age);
	  double Wsnew = rho*Asnew*GetValue(*ts,LGAL);
	  SetValue(data,LGAAs,Asnew+GetValue(*ts,LGAAs));
	  //Mass of sapwood used in diamater growth
	  SetValue(data,LGAiWs,GetValue(data,LGAiWs)+Wsnew);
	  //Total foliage
	  SetValue(data,DGWf,GetValue(data,DGWf)+GetValue(*ts,LGAWf));
	}
      }
      return data;
    }
    else {//(mode == LGMGROWTH)
      //Dimensions of the new segments (age == 0) are iteratively set.
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	if (GetValue(*ts,LGAage) == 0.0){
	  //cout << "Diameter Growth New Segment " << GetValue(*ts, LGAomega) << " "  << GetValue(data,LGAomega) <<endl;
	  //Mass of the new sapwood
	  double Wsnew = GetValue(*ts,LGAWs);
	  SetValue(data,LGAAs,GetValue(*ts,LGAAs));
	  //Collect masses for bookkeeping
	  SetValue(data,LGAiWf,GetValue(*ts,LGAWf));
	  SetValue(data,DGWf,GetValue(*ts,LGAWf));
	  SetValue(data,LGAiWs,Wsnew);
	}
	//Set the dimensions of the existing segments
	else{// (GetValue(*ts,LGAage) > 0.0)
	  //cout << "Diameter Growth Old Segment " << GetValue(*ts, LGAomega) 
	  //<< " "  << GetValue(data,LGAomega) <<endl;
	  const ParametricCurve& fm = GetFunction(GetTree(*ts),LGMFM);
	  //Sapwood requirement of  remaining foliage, assume fm returns
	  //proportion initial foliage  present, declining function from
	  //1 to 0.
	  LGMdouble Asr = fm(GetValue(*ts,LGAage))*GetValue(*ts,LGAAs0);
	  LGMdouble Asu = GetValue(data,LGAAs); //sapwood area from above
	  LGMdouble Ahown  = GetValue(*ts,LGAAh);//own heartwood
	  //possible new radius
	  LGMdouble Rnew = sqrt((Asu + Ahown + Asr)/PI_VALUE);
	  //compare Rnew to R, choose max
	  Rnew = max(Rnew, GetValue(*ts,LGAR));
	  //New sapwood requirement, thickness growth
	  double Asnew = PI_VALUE*pow(Rnew,2.0) -  GetValue(*ts,LGAA);
	  //New wood radius
	  SetValue(*ts,LGAR,Rnew);
	  //Mass of  the new  sapwood This is  a shorthand  fix: Hakkila
	  //1969 CIFF 67.6 shows that wood density for branches is > 200
	  //kgC/m3. To achieve this add 10 years to age, if branch
	  const ParametricCurve& wd = 
	    GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*tc)), SPWD);
	  double age = GetValue(*tc, LGAage);
	  if (GetValue(*ts,LGAomega) > 1.0)
	    age += 10.0;	
	  double rho = wd(age);
	  double Wsnew = rho*Asnew*GetValue(*ts,LGAL);
	  SetValue(data,LGAAs,GetValue(*ts,LGAAs));
	  //Mass of sapwood used in diamater growth
	  SetValue(data,LGAiWs,GetValue(data,LGAiWs)+Wsnew);
	  //Total foliage
	  SetValue(data,DGWf,GetValue(data,DGWf)+GetValue(*ts,LGAWf));
	}
      }
      return data;
    }
  }
private:
  LGMALLOCATORMODE mode;
};
#endif
