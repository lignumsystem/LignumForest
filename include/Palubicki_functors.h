#ifndef PALUBICKI_FUNCTORS_H
#define PALUBICKI_FUNCTORS_H
#include <Lignum.h>
///\file Palubicki_functors.h
///\brief Functors accoring to Self-organizing tree models for image synthesis
///
///W. Palubicki and K. Horel and S. Longay and
///A. Runions and B. Lane and R. Mech and P. Prusinkiewicz. 2009.
///*Self-organizing tree models for image synthesis* ACM Transactions on
///Graphics 28 58:1-10.
namespace LignumForest{
  ///This is according to W. Palubicki and K. Horel and S. Longay and
  ///A. Runions and B. Lane and R. Mech and P. Prusinkiewicz. 2009.
  ///Self-organizing tree models for image synthesis ACM Transactions on
  ///Graphics 28 58:1-10. Length growth is according to
  ///Extended Borchert-Honda model, BH_resource is share [0,1] of total
  ///intercepted radiation that has been is here distributed acropetally 
  ///according to their function.
  class EBH_basipetal_info {
  public:
    EBH_basipetal_info() : orders(1,1.0), Qvs(1,0.0) {;}
    EBH_basipetal_info& operator += (const EBH_basipetal_info& add_this) {
      for(int i = 0; i < (int)(add_this.orders).size(); i++) {
	orders.push_back((add_this.orders)[i]);
      }
      for(int i = 0; i < (int)(add_this.Qvs).size(); i++) { //could be (add_this.orders).size()
	Qvs.push_back((add_this.Qvs)[i]);              //they are the same
      }
      return *this;
    }
    vector<double> orders;
    vector<double> Qvs;
  };

  class EBH_basipetal {
  public:
    EBH_basipetal(const ParametricCurve la_in, const int mode) : lambda_fun(la_in), ebh_mode(mode) {;}
    EBH_basipetal_info& operator () (EBH_basipetal_info& val, 
				     TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const {

      if(ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)) {
	(val.orders)[0] = GetValue(*ts, LGAomega); //vectors orders, Qvs are of length 1 here
	//      (val.Qvs)[0] += GetValue(*ts, LGAQin);
	switch(ebh_mode) {
	case 1:
	  (val.Qvs)[0] += GetValue(*ts, LGAQabs);
	  break;
	case 2:
	  (val.Qvs)[0] += GetValue(*ts, SPrue) * GetValue(*ts, LGAQabs);
	  break;
	default:
	  (val.Qvs)[0] += GetValue(*ts, LGAQin);
	}
      
	ts->setQv((val.Qvs)[0]);
      }

      if(BranchingPoint<ScotsPineSegment,ScotsPineBud>* bp = 
	 dynamic_cast<BranchingPoint<ScotsPineSegment,ScotsPineBud>*>(tc)) {

	double sum = 0.0;
	double my_order = GetValue(*bp,LGAomega);
	for(int i = 0; i < (int)(val.orders).size(); i++) {
	  if((val.orders)[i] > my_order) {
	    sum += (1.0 - lambda_fun(my_order)) * (val.Qvs)[i];
	  } else {
	    sum += lambda_fun(my_order) * (val.Qvs)[i];
	  }
	}
	SetValue(*bp, LGAMaxD,sum);

	sum = accumulate((val.Qvs).begin(),(val.Qvs).end(),0.0);
	(val.orders).clear();
	(val.Qvs).clear();
	(val.orders).push_back(my_order);
	(val.Qvs).push_back(sum);
      }

      return val;
    }

  private:
    ParametricCurve lambda_fun;
    int ebh_mode;
  };


  class EBH_acropetal_info  {
  public:
    EBH_acropetal_info(const double o_in, const double r_in, const double s_in) : 
      order(o_in), resource(r_in), sum(s_in) {;}
    double order;
    double resource;
    double sum;
  };



  class EBH_acropetal{
  public:
    EBH_acropetal(const ParametricCurve la_in) : lambda_fun(la_in) {}
    EBH_acropetal_info& operator()(EBH_acropetal_info& coming_up, 
				   TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
    {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	double my_order = GetValue(*ts, LGAomega);
	double my_resource = 0.0;
	//      cout << "TS order point " << my_order << " " << GetPoint(*ts);

	//The cumulative Qin stored at BP below this TS may be == 0, if this
	//is a TS at the end of an axis with R == 0 & L == 0 (a TS that 
	//had in the iteration L < Lmin and therefore was set L = 0).
	//coming_up.sum = cumulative Qin stored at BP below this TS.
	//Cannot divide by 0, set my_resource = 0.
	double cu_resource = coming_up.sum;
	if(cu_resource < R_EPSILON) {
	  my_resource = 0.0;
	} else {
	  if(my_order > coming_up.order){
	    my_resource = coming_up.resource * (1.0 - lambda_fun(coming_up.order))
	      * ts->getQv() / coming_up.sum;  //lateral
	  } else{
	    my_resource = coming_up.resource * lambda_fun(coming_up.order)
	      * ts->getQv() / coming_up.sum;
	  }
	}
	ts->setEBHResource(my_resource);

	coming_up.resource = my_resource;
      }
      if(BranchingPoint<ScotsPineSegment,ScotsPineBud>* bp = 
	 dynamic_cast<BranchingPoint<ScotsPineSegment,ScotsPineBud>*>(tc)) {
	coming_up.sum = GetValue(*bp,LGAMaxD);
	coming_up.order = GetValue(*bp, LGAomega);
      }
      return coming_up;
    }

  private:
    ParametricCurve lambda_fun;
  };


  class MaxEBHResource_info {
  public:
    MaxEBHResource_info(): my_resource(-R_HUGE) {}
    MaxEBHResource_info& operator += (const MaxEBHResource_info& add_this) {
      if(add_this.my_resource > my_resource) {
	my_resource = add_this.my_resource;
      }
      return *this;
    }
    double my_resource;
  };

  class MaxEBHResource{
  public:
    MaxEBHResource_info& operator()(MaxEBHResource_info& value, 
				    TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const {
      if (Axis<ScotsPineSegment,ScotsPineBud>* ax =
	  dynamic_cast<Axis<ScotsPineSegment,ScotsPineBud>*>(tc)){
	TreeSegment<ScotsPineSegment,ScotsPineBud>* ts = GetLastTreeSegment(*ax);

	if(ts != NULL) {
	  value.my_resource = max((dynamic_cast<ScotsPineSegment*>(ts))->getEBHResource(),
				  value.my_resource);
	}
      }
      return value;
    }
  };

  class NormalizeEBHResource{
  public:
    NormalizeEBHResource(const double mv) : max_value(mv) {}
    TreeCompartment<ScotsPineSegment,ScotsPineBud>*
    operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
    {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	ts->setEBHResource(ts->getEBHResource()/max_value);
      }

      return tc;
    }
  private:
    double max_value;
  };


  class printEBHR {
  public:
    TreeCompartment<ScotsPineSegment,ScotsPineBud>* operator () ( 
								 TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const {
      if (Axis<ScotsPineSegment,ScotsPineBud>* ax =
	  dynamic_cast<Axis<ScotsPineSegment,ScotsPineBud>*>(tc)){
	TreeSegment<ScotsPineSegment,ScotsPineBud>* ts = GetLastTreeSegment(*ax);

	if(ts != NULL) {
	  cout << "E " << (dynamic_cast<ScotsPineSegment*>(ts))->getEBHResource() << " " << GetPoint(*ts);
	}
      }
      return tc;
    }
  };

  class printqin {
  public:
    TreeCompartment<ScotsPineSegment,ScotsPineBud>* operator () ( 
								 TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)) {

	cout << "Q " << GetValue(*ts, LGAQin) << endl;
      }
      return tc;
    }
  };
}//End namespace LignumForest
#endif
