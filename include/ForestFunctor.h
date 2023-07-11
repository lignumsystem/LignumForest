#ifndef FOREST_FUNCTORS_H
#define FOREST_FUNCTORS_H

#include <ScotsPine.h>
#include <SomeFunctors.h>
#include <sstream>
using namespace std;

namespace LignumForest{
  ///Collect   new   foliage   and   set   initial   root   mass.    Move
  ///CollectNewFoliage to stl-lignum.
  class InitialRootMass{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      double wf = 0.0;
      Accumulate(*t,wf,CollectNewFoliage<ScotsPineSegment,ScotsPineBud>());
      SetValue(*t,TreeWr,0.0);
      SetValue(*t,TreeWr,GetValue(*t,LGPar)*wf);
    }
  };

  ///Set STARmean  to segments 
  class ForestSetStar{
  public:
    ForestSetStar(const SetStarMean<ScotsPineSegment,ScotsPineBud>& sm)
      :starmean(sm){}
    void  operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      ForEach(*t,starmean);
    }
  private:
    const SetStarMean<ScotsPineSegment,ScotsPineBud>& starmean;
  };

  class ForestDumpFoliage{
  public:
    ForestDumpFoliage(VoxelSpace& voxel_space,int parts)
      :vs(voxel_space),num_parts(parts){}
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      try{
	DumpScotsPineTree(vs,*t,num_parts);
      }
      catch (OutOfVoxelSpaceException e){
	ostringstream error;
	Point box = e.getBox();
	Point p = e.getPoint();
	Point tp = GetPoint(*t);
	error << "Foliage dump exception with point " 
	      << p.getX() << " " << p.getY() << " " << p.getZ() << endl;
	error << "Box indices " << box.getX() << " " << box.getY() << " "
	      << box.getZ() << endl;
	error << "Tree Point " << tp.getX() << " " << tp.getY() << " " 
	      << tp.getZ() <<endl;
	LGMError(error.str());
	//Remove the tree from  simulation. We cannot simply say: delete
	//t; t=NULL;  it does  not work. We  instead mark the  tree dead
	//(age < 0)  and remove the tree, its  L-system and its position
	//in RemoveDeadTrees().
	SetValue(*t,LGAage,-100.0);
      }
    }
  private:
    VoxelSpace& vs;
    int num_parts;
    int n;
  };

  class ForestPairwiseLight{
  public:
    ForestPairwiseLight(EvaluateRadiationForCfTreeSegment<ScotsPineSegment,ScotsPineBud>& Rad)
      :R(Rad){}
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      ForEach(*t,R);
    }
  private:
    EvaluateRadiationForCfTreeSegment<ScotsPineSegment,ScotsPineBud>& R;
  };

  class SetForestTreeQabs{
  public:
    SetForestTreeQabs(VoxelSpace& voxel_space,int parts)
      :vs(voxel_space),num_parts(parts){}
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      SetScotsPineTreeQabs(vs,*t,num_parts); 
    }
    VoxelSpace& vs;
    int num_parts;
  };

  class ForwardForestWf{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      double wf = 0.0;
      PropagateUp(*t,wf,ForwardWf<ScotsPineSegment,ScotsPineBud>());
    }
  };
    
  class ForestPhotosynthesis{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      t->photosynthesis();
    }
  };

  class ForestRespiration{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      t->respiration();
    }
  };

  class ForestAging{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      //Aging (see actual implementation for segment method aging())
      ForEach(*t,TreeAging<ScotsPineSegment,ScotsPineBud>());
      //Root mortality
      SetValue(*t,TreeWr, 
	       GetValue(*t,TreeWr)-GetValue(*t,LGPsr)*GetValue(*t,TreeWr));
    }
  };

  class ForestForwardQin{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      double qin = 0.0;
      PropagateUp(*t,qin,ForwardQin<ScotsPineSegment,ScotsPineBud>());
    }
  };

  class ForestVigourIndex{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      TreePhysiologyVigourIndex(*t);
    }
  };

  //QinMax relative to sky (ball sensor reading)
  class ForestPairwiseQinMax{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      SetValue(*t,TreeQinMax,GetFirmament(*t).diffuseBallSensor());
    }
  };

  //QinMax relative to tree itself
  class ForestVoxelSpaceQinMax{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      LGMdouble qin_max = 0.0;
      qin_max = Accumulate(*t,qin_max,GetQinMax<ScotsPineSegment,ScotsPineBud>());
      SetValue(*t,TreeQinMax,qin_max);
    }
  };

  class ForestPathLength{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      LGMdouble plength = 0.0;
      PropagateUp(*t,plength,PathLength<ScotsPineSegment,ScotsPineBud>());
    }
  };

  class ForestKillBuds{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      double qin=0.0;
      PropagateUp(*t,qin,KillBuds<ScotsPineSegment,ScotsPineBud>());
    }
  };

  class ForestAllocateGrowth{
  public:
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      LPineGrowthFunction<ScotsPineSegment,ScotsPineBud> G(*t);
      try{
	G.init();
	Bisection(0.0,5.0,G,0.01);
      }//The bracketing did not succeed
      catch(BisectionBracketException be){
	ostringstream error;
	Point p = GetPoint(*t);
	error << "BracketException " << be.getFa() << " " << be.getFb() << " "
	      << be.getFbl() << " "  << be.getA() << " " << be.getB() << " " 
	      << be.getBl() << endl;
	error << "Tree Point " << p.getX() << " " << p.getY() 
	      << " " << p.getZ() << endl;
	LGMError(error.str());
	//Remove the tree from  simulation. We cannot simply say: delete
	//t; t=NULL;  it does  not work. We  instead mark the  tree dead
	//(age < 0)  and remove the tree, its  L-system and its position
	//in RemoveDeadTrees().
	SetValue(*t,LGAage,-100.0);
      }//The root finding did not succeed
      catch(BisectionMaxIterationException mie){
	ostringstream error;
	Point p = GetPoint(*t);
	error << "MaxIterationException " <<mie.getFa()<<" "<<mie.getFb()<<" "
	      << mie.getFc() << " "  << mie.getA() << " " << mie.getB() << " " 
	      << mie.getC() << endl;
	error << "tree Point " << p.getX() << " " << p.getY() 
	      << " " << p.getZ() << endl;
	LGMError(error.str());
	//Remove the tree from  simulation. We cannot simply say: delete
	//t; t=NULL;  it does  not work. We  instead mark the  tree dead
	//(age < 0)  and remove the tree, its  L-system and its position
	//in RemoveDeadTrees().
	SetValue(*t,LGAage,-100.0);
      }//P < M --> nothing to allocate
      catch(TreeGrowthAllocatorException ae){
	ostringstream error;
	Point p = GetPoint(*t);
	error << "P < M " << ae.getP() << " " << ae.getM() << endl;
	error << "Tree Point: " <<  p.getX() << " " << p.getY() 
	      << " " << p.getZ() <<endl;
	LGMError(error.str());
	//Remove the tree from  simulation. We cannot simply say: delete
	//t; t=NULL;  it does  not work. We  instead mark the  tree dead
	//(age < 0)  and remove the tree, its  L-system and its position
	//in RemoveDeadTrees().
	SetValue(*t,LGAage,-100.0);
      }
      //If the tree is still alive
      if (GetValue(*t,LGAage) >= 0){
	//Kill the buds (LGAL < R_EPSILON)
	bool kill = false;
	PropagateUp(*t,kill,KillBudsAfterAllocation<ScotsPineSegment,ScotsPineBud>());
	//Bisection (allocation) done, now do the diameter growth
	TcData data;
	AccumulateDown(*t,data,TreeDiameterGrowth<ScotsPineSegment,ScotsPineBud>());
	double wfnew = 0.0;
	Accumulate(*t,wfnew,CollectNewFoliage<ScotsPineSegment,ScotsPineBud>());
	//Update root mass 
	SetValue(*t,TreeWr, 
		 GetValue(*t,TreeWr)+
		 GetValue(*t,LGPar)*wfnew);
      }
    }
  };

  class ForestDiffuseRadiation{
  public:
    ForestDiffuseRadiation(ParametricCurve& K1):K(K1){}
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
      //Diameter and heigth at the crown base.
      DCLData dcl;
      AccumulateDown(*t,dcl,AddBranchWf(),DiameterCrownBase<ScotsPineSegment,ScotsPineBud>());
      //Collect foliage 
      LGMdouble wf = 0.0;
      Accumulate(*t,wf,CollectFoliageMass<ScotsPineSegment,ScotsPineBud>());
      LGMdouble Hc = dcl.HCrownBase();
      LGMdouble Af = GetValue(*t,LGPsf)*wf;
      DiffuseForestRadiation<ScotsPineSegment,ScotsPineBud> FRad(Hc,Af,K);
      ForEach(*t,FRad);
    }
  private:
    ParametricCurve& K;
  };

  //Print for  each tree some summary information:  length, trunk shape,
  //crown shape
  class ForestPrintSummary{
  public:
    ForestPrintSummary(ofstream& file);
    void operator()(Tree<ScotsPineSegment,ScotsPineBud>* t);
  private:
    ofstream& f;
  };

  //Print to  file for each tree  1)Position, 2) Length,  3) Dbase, 4)
  //Dbh, 5)  Height crown  base, 6) Diameter  crown base and  7) Crown
  //height
  inline ForestPrintSummary::ForestPrintSummary(ofstream& file)
    :f(file)
  {
    f << setfill(' ');
    f << setw(11) << "X" << setw(11) << "Y" << setw(11) << "H" 
      << setw(11) << "Dbase" << setw(11) << "Dbh" << setw(11) << "Hcb" 
      << setw(11) << "Dcb" << setw(11) << "Hc" << endl;
  }

  inline void ForestPrintSummary::operator()(Tree<ScotsPineSegment,ScotsPineBud>* t){
    Point p = GetPoint(*t);
    LGMdouble h = GetValue(*t,LGAH);
    LGMdouble dbase = GetValue(*t,LGADbase);
    LGMdouble dbh = GetValue(*t,LGADbh);
    DCLData dcl;
    AccumulateDown(*t,dcl,AddBranchWf(),DiameterCrownBase<ScotsPineSegment,ScotsPineBud>());
    LGMdouble hcb = dcl.HCrownBase();//height crown base
    LGMdouble dcb = dcl.DCrownBase();//diameter crown base
    LGMdouble hc = h-hcb; //Crown height
    
    f << setfill(' ');
    f << setw(11) << p.getX() << setw(11) << p.getY() 
      << setw(11) << h << setw(11) << dbase << setw(11) << dbh 
      << setw(11) << hcb << setw(11) << dcb << setw(11) << hc << endl; 
  }
}//End namespace LignumForest
#endif
