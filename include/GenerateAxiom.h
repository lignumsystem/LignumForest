#ifndef GENERATE_AXIOM_H
#define GENERATE_AXIOM_H
#include <Lignum.h>
#include <fstream>
using namespace std;

namespace LignumForest{
  ///Track SB() and EB() balance
  static int level = 0;
  ///Generate the beginning of the axiom. The beginning is designed  for the
  ///ScotsPine used in LignumForest
  void GenerateBeginAxiom(fstream& file)
  {
    file << "Start:"<< endl << "{" << endl 
	 << "PineBudData data(ALIVE,0.0,1.0,0.1);" << endl
	 << "mode = 0;" << endl <<endl;
    file << "produce " << flush;
  }

  ///End of Axiom adds the semicolon and right curly brace
  void GenerateEndAxiom(fstream& file)
  {
    file << ";" << endl << "}" << flush;
  }

  template <class TS,class B>
  void GenerateAxiom(fstream& file, Tree<TS,B>& t)
  {
    GenerateBeginAxiom(file);
    GenerateAxiom(file,&GetAxis(t));
    GenerateEndAxiom(file);
  }

  ///Generate  axiom. The  idea  is to  move  to Turtle  with MoveTo  and
  ///SetHead commands. This is how we do not have to calculate rotations.
  ///\note Declare  MoveTo(double,double,double) and
  ///SetHead(double,double,double) modules manually in the L-system.
  template <class TS,class B>
  void GenerateAxiom(fstream& file, TreeCompartment<TS,B>* tc)
  {
    if (TS* ts = dynamic_cast<TS*>(tc)){
      Point p = GetPoint(*ts);
      Point e = GetDirection(*ts);
      file << "MoveTo("<<p.getX()<<","<<p.getY()<<","<<p.getZ()<<")"<< flush;
      file << "SetHead("<<e.getX()<<","<<e.getY()<<","<<e.getZ()<<")"<<flush;
      file << "F("<<GetValue(*ts,LGAL)<<")"<<flush;
    }
    else if (BranchingPoint<TS,B>* bp = dynamic_cast<BranchingPoint<TS,B>*>(tc)){    
      list<Axis<TS,B>*>& axis_ls = GetAxisList(*bp);
      typename list<Axis<TS,B>*>::iterator first = axis_ls.begin();
      typename list<Axis<TS,B>*>::iterator last = axis_ls.end();
      if (axis_ls.empty())
	file << "SB()EB()" << flush;
      while (first != last){
	file << "SB()" << flush;
	level++;
	GenerateAxiom(file,*first);
	file << "EB()" <<flush;
	level--;
	//This is for the ease of reading the possibly long axiom file
	if (level == 0)
	  file << endl << endl;
	first++;
      }

    }
    else if (B* b = dynamic_cast<B*>(tc)){
      Point p = GetPoint(*b);
      Point e = GetDirection(*b);
      double omega = GetValue(*b,LGAomega);
      file << "MoveTo("<<p.getX()<<","<<p.getY()<<","<<p.getZ()<<")"<< flush;
      file << "SetHead("<<e.getX()<<","<<e.getY()<<","<<e.getZ()<<")"<<flush;
      file << "B(data,"<<omega<<",1.0)" << flush;
    }
    else if (Axis<TS,B>* axis = dynamic_cast<Axis<TS,B>*>(tc)){
      list<TreeCompartment<TS,B>*>& tc_ls = GetTreeCompartmentList(*axis);
      typename list<TreeCompartment<TS,B>*>::iterator first = tc_ls.begin();
      typename list<TreeCompartment<TS,B>*>::iterator last = tc_ls.end();
      while (first != last){
	GenerateAxiom(file,*first++);
      }
    }
    else{
    }
  }
}//End namespace LignumForest
#endif
