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
	  ///\par Branch length
	  ///Branch length as a sum of segment lengths
	  ///\internal
	  ///\snippet{lineno} branchfunctor.cc LB
	  // [LB]
	  while(I != sl.end()) {
	    if (ScotsPineSegment* seg = 
		dynamic_cast<ScotsPineSegment*>(*I)) {
	      lb += GetValue(*seg, LGAL);
	    }
	    I++;
	  }
	  // [LB]
	  ///\endinternal
	  ///\par Diameter data
	  ///\internal
	  ///\snippet{lineno} branchfunctor.cc DD
	  //[DD]
	  //Diameter of the first segment 
	  double d = 2.0 * GetValue(*fs, LGAR);
	  id.d2 += d * d;
	  id.d2l += d * d * lb;
	  id.lsum += lb;
	  id.n_br++;
	  //[DD]
	  ///\endinternal
	}
      }//if fs != NULL
    }//if (Axis....)
    return id;
  }
}//End namespace LignumForest
