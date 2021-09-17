#ifndef STANDDESCRIPTORI_H
#define STANDDESCRIPTORI_H
#include <StandDescriptor.h>

template<class TREE>
void StandDescriptor<TREE>::evaluateStandVariables(vector<TREE*> vt, vector<pair<double,double> >& loc)
{
  dbh_mean = 0.0; hdom = -R_HUGE; stemVolume = 0.0;
  standBasalArea = 0.0; meanHeight = 0.0; basalAreaAtCrownBase = 0.0; minDBH = R_HUGE;
  maxDBH = -R_HUGE; minH = R_HUGE; maxH = -R_HUGE; 
  age = -INT_MAX; noTrees = 0.0;
  dbase_mean = 0.0; minDbase = R_HUGE; maxDbase = -R_HUGE;
  

  double d = 0, d2S_dbh = 0.0, h = 0.0, d2S_base = 0.0;
  double dsum_base = 0.0, dsum_dbh = 0.0;
  double CL_sum = 0.0;
  
  double d2_dbh, d2_base;
  stemVolume = standBasalArea = 0.0;
  basalAreaAtCrownBase = 0.0;
  minCrownLimit = R_HUGE; 	
  maxCrownLimit = -R_HUGE;
  meanCrownLimit = 0.0;
  LAI = 0.0;
  Wf = 0.0;

  int vector_size = (int)vt.size();
  noTrees = 0;
   
  for(unsigned int i = 0; i < (unsigned int)vector_size; i++) {
    Point p(loc[i].first, loc[i].second,0.0);
    if(inPlot(p)){
      noTrees++;
      int a = (int)GetValue(*vt[i],LGAage);  //
      if(a > age) age = a;
      d = GetValue(*vt[i],LGADbh);
      if(d < minDBH) minDBH = d;
      if(d > maxDBH) maxDBH = d;
      d2_dbh = d * d;
      d2S_dbh += d2_dbh;
      dsum_dbh += d * d2_dbh;
      standBasalArea += PI_VALUE * d2_dbh / 4.0;

      d = GetValue(*vt[i],LGADbase);
      if(d < minDbase) minDbase = d;
      if(d > maxDbase) maxDbase = d;
      d2_base = d * d;
      d2S_base += d2_base;
      dsum_base += d * d2_base;

      h = GetValue(*vt[i],LGAH);
      if(h > hdom) hdom = h;
      if(h > maxH) maxH = h;
      if(h < minH) minH = h;
      meanHeight += h * d2_base;                     //OK, painotetaan base diameterilla

      stemVolume +=  GetValue(*vt[i],LGAWstem)/ GetValue(*vt[i],LGPrhoW);  //Tama korjattava!

      DCLData dcl;
      dcl.DCrownBase(d);   //Set values of base of tree to get it start correctly (d = base diameter)
      dcl.HCrownBase(0.0);

      AccumulateDown(*vt[i],dcl,AddBranchWf(),DiameterCrownBase<ScotsPineSegment,ScotsPineBud>());
      LGMdouble hcb = dcl.HCrownBase();
      if(hcb < minCrownLimit)	minCrownLimit = hcb;
      if(hcb > maxCrownLimit)	maxCrownLimit = hcb;
      CL_sum += d2_dbh * hcb;
      

      LGMdouble dcb = dcl.DCrownBase();
      basalAreaAtCrownBase += dcb*dcb;

      LGMdouble wf = 0.0;
      Wf += Accumulate(*vt[i], wf, CollectFoliageMass<ScotsPineSegment,ScotsPineBud>());

      LGMdouble af = 0.0;
      LAI += Accumulate(*vt[i], af, CollectFoliageArea<ScotsPineSegment,ScotsPineBud>());
    }  //if(this->inPlot(...
  } //for(unsigned int i = no_trees ...

			
  if(d2S_dbh > 0.0)	{
    dbh_mean = dsum_dbh / d2S_dbh;
    meanCrownLimit =  CL_sum / d2S_dbh;
  }
  else	{
    dbh_mean = 0.0;
    meanCrownLimit = 0.0;
  }
  if(d2S_base > 0.0)	{
    dbase_mean = dsum_base / d2S_base;
    meanHeight /= d2S_base;
  }
  else	{
   dbase_mean = 0.0;
   meanHeight = 0.0;
  }

  if(area > 0.0) {
  stemVolume /= area;
  standBasalArea /= area;
  basalAreaAtCrownBase *= PI_VALUE/(4.0*area);
  LAI /= area;
  Wf /= area;
  }

}


template<class TREE>
void StandDescriptor<TREE>::writeOutput(ofstream& stand_output) {

  stand_output << 10000.0*(double)noTrees/area << " " << dbase_mean << " " << minDbase << " "
	       << maxDbase << " " << dbh_mean << " "
	       << minDBH << " " << maxDBH << " " << meanHeight << " " << minH << " " << maxH << " "
	       << standBasalArea << " " << basalAreaAtCrownBase << " " << stemVolume << " "
	       << LAI << " " << Wf << " " << meanCrownLimit << endl;
}



#endif
