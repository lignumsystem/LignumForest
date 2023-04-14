///\file  branchfunctor.cc
///\brief Collect simple descriptive data from branches.
///
///Collect descriptive data from branches with Gravelius order == 2.
#include <SomeFunctors.h>
namespace LignumForest{
///Branch data for Gravelius order 2
summing& Branchmeans::operator()(summing& id,
				 TreeCompartment<ScotsPineSegment, ScotsPineBud>*
				 tc)const {
  TreeSegment<ScotsPineSegment,ScotsPineBud>* fs = NULL;
  if (Axis<ScotsPineSegment, ScotsPineBud>* axis =
      dynamic_cast<Axis<ScotsPineSegment, ScotsPineBud>*>(tc)){
    fs = GetFirstTreeSegment(*axis);
    if (fs != NULL){
      if(GetValue(*fs, LGAomega) == 2) {
	std::list<TreeCompartment<ScotsPineSegment, ScotsPineBud>*>&
	  sl = GetTreeCompartmentList(*axis);
	double lb = 0.0;
	list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>::iterator I
	  = sl.begin();
	///\section branchmeans Collecting branch data
	///\subsection lbdata Branch length
	///Branch length as a sum of segment lengths
	///\snippet{lineno} branchfunctor.cc LB
	///\internal
	//[LB]
	while(I != sl.end()) {
	  if (ScotsPineSegment* seg = 
	      dynamic_cast<ScotsPineSegment*>(*I)) {
	    lb += GetValue(*seg, LGAL);
	  }
	  I++;
	}
	//[LB]
	///\endinternal
	///\subsection diam Diameter data
	///Collect branch data
	///\snippet{lineno} branchfunctor.cc DD
	///\internal
	//[DD]
	//Diameter of the first segment 
	double d = 2.0 * GetValue(*fs, LGAR);
	//Sum of diameter squared from the first segments ("area")
	id.d2 += d * d;
	//Sum of diameter squared times length ("volume")
	id.d2l += d * d * lb;
	//Sum of length of branches
	id.lsum += lb;
	//Sum of branches Gravelius order 2
	id.n_br++;
	//[DD]
	///\endinternal
      }
    }//if fs != NULL
  }//if (Axis....)
  return id;
}
}
