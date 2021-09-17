#include <SomeFunctors.h>
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
	while(I != sl.end()) {
	  if (ScotsPineSegment* seg = 
	      dynamic_cast<ScotsPineSegment*>(*I)) {
	    lb += GetValue(*seg, LGAL);
	  }
	  I++;
	}
	double d = 2.0 * GetValue(*fs, LGAR);
	id.d2 += d * d;
	id.d2l += d * d * lb;
	id.lsum += lb;
	id.n_br++;
      }
    }//if fs != NULL
  }//if (Axis....)
  return id;
}
