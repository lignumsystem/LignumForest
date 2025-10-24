#ifndef SENSITIVITYI_H
#define SENSITIVITYI_H
#include <Sensitivity.h>
/// \file SensitivityI.h
/// \brief Print tree data for sensitivity analysis impelmentation
/// \deprecated Use HDF5 files for data analysis
namespace LignumForest{
template <class TS,class BUD>
void Sensitivity<TS,BUD>::printHeader(const string& file_name)
{
  file.open(file_name.c_str(),ofstream::app);
  file << "Age       Dbase      H      H/Dbase       CL        DCL       Hc      " << flush;
  file << "Abase        Aswbase         Ahwbase         " <<flush;
  file << "Aclimit      Aswclimit          Ahwclimit             " <<flush;
  file << "Wood        Wsw        Whw      Wstem       Wbranch        " << flush;
  file << "Wf          Wf0        Wr" << endl;
}

template <class TS,class BUD>
Sensitivity<TS,BUD>::~Sensitivity()
{
  file.close();
}

template <class TS,class BUD>
void Sensitivity<TS,BUD>::printData(Tree<TS,BUD>& t)
{
  //diameter and height at crown base
  DCLData dcl;
  AccumulateDown(t,dcl,AddBranchWf(),DiameterCrownBase<TS,BUD>());

  //Collect all foliage 
  LGMdouble wfoliage = 0.0;
  Accumulate(t,wfoliage,CollectFoliageMass<TS,BUD>());

  //Collect new foliage (segment age = 0)
  double wf0 = 0.0;
  Accumulate(t,wf0,CollectNewFoliageMass<TS,BUD>());

  //Collect wood mass in a tree
  double wood = 0.0;
  Accumulate(t,wood,CollectWoodMass<TS,BUD>());

  //Collect wood mass in the main axis
  double wstem = GetValue(t,LGAWstem);

  //Wood mass in branches
  double wbranches = wood - wstem;
  
  //Collect sapwood mass (crown + stem)
  double ws =0.0;
  Accumulate (t,ws,CollectSapwoodMass<TS,BUD>());

  double wh = 0.0;
  //Collect heartwood mass
  Accumulate (t,wh,CollectHeartwoodMass<TS,BUD>());

  LGMdouble age = GetValue(t,LGAage);    //Tree age 
  LGMdouble dbase = GetValue(t,LGADbase);//diameter base
  LGMdouble h = GetValue(t,LGAH);        //tree height
  LGMdouble h_dbase_ratio = h/dbase;    //height / diameter_base
  LGMdouble cl  = dcl.HCrownBase();      //crown base height
  LGMdouble cld = dcl.DCrownBase();      //diameter crown base
  LGMdouble hc = h - cl;                 //crown height
  LGMdouble abase = GetValue(t,LGAAbase);//area at base
  LGMdouble aswbase = GetValue(t,LGAAsbase);//sapwood area at base
  LGMdouble ahwbase = GetValue(t,LGAAhwbase);//heartwood area at base
  LGMdouble aclimit = dcl.ACrownBase();      //area at crown base
  LGMdouble aswclimit = dcl.ASwCrownBase();  //sapwood area at crown base
  LGMdouble ahwclimit = dcl.AHwCrownBase();  //heartwood area at crown base
  LGMdouble wr = GetValue(t,TreeWr);      //root mass

  file << left << setfill(' ')
       << setw(14) << age
       << setw(14) << dbase 
       << setw(14) << h
       << setw(14) << h_dbase_ratio
       << setw(14) << cl
       << setw(14) << cld
       << setw(14) << hc
       << setw(14) << abase
       << setw(14) << aswbase
       << setw(14) << ahwbase
       << setw(14) << aclimit
       << setw(14) << aswclimit
       << setw(14) << ahwclimit
       << setw(14) << wood
       << setw(14) << ws
       << setw(14) << wh
       << setw(14) << wstem
       << setw(14) << wbranches
       << setw(14) << wfoliage
       << setw(14) << wf0
       << setw(14) << wr <<endl;
}
}//End namespace LignumForest
#endif
