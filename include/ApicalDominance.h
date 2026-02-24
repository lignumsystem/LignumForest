#ifndef APICAL_DOMINANCE_H
#define APICAL_DOMINANCE_H
#include <Uniform.h>
///\file ApicalDominance.h
///\brief Experimental functors to control apical dominance in growth
///
///During later stage in development Scots pine tend to lose the strict apical
///dominance of younger trees.

namespace LignumForest{
  ///\brief Supress growth in the top part of the crown
  ///\tparam TS Tree segment
  ///\tparam BUD Bud
  template <class TS,class BUD>
  class SupressCrownTop{
  public:
    ///\brief Constructor
    ///\param t_height Tree height defined as the highest point in the tree
    ///\param r_height Relative height in the interval [0,1] defining the top part of the crown
    ///\param p Probability of a bud in the top part of the crown to die
    ///\param seed Seed for the pseudo random number generator \sa SupressCrownTop::u
    ///\pre \f$ 0 \leq \f$ \p r_height \f$ \leq 1 \f$
    ///\pre \p seed < 0 initializes the sequence of pseudo random numbers \sa cxxadt::Uniform
    ///\post The sequence of pseudo random numbers initialized 
    SupressCrownTop(double t_height,double r_height,double p,int seed=-1)
      :tree_height(t_height),crown_top(r_height),p_dead(p){u(seed);}
    ///\brief Copy constructor
    ///\param sct SupressCrownTop instance
    ///\pre SupressCrownTop::u in \p sct has been initialized.
    SupressCrownTop(const SupressCrownTop& sct)
      :tree_height(sct.tree_height),crown_top(sct.crown_top),p_dead(sct.p_dead){}
    ///\brief Assignment
    ///\param sct SupressCrownTop instance
    ///\return self  
    SupressCrownTop& operator=(const SupressCrownTop& sct);
    ///\brief Supress tree growth
    ///
    ///Set a bud DEAD with a probability SupressCrownTop::p_dead in the top part of the crown
    ///\param tc Tree compartment
    ///\retval tc The tree compartment 
    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;
  private:
    double tree_height;///< Highest z-coordinate in the tree
    ///Top part in the interval [0,1] of the tree crown.
    ///For example buds \f$ 0.9 \times \mathrm{tree_height}\f$ and above are in the top part of the crown.
    double crown_top;
    double p_dead;///< Probability of a bud in the top part of the crown to die
    ///Uniformly distributed pseudo random number generator
    ///\note The  random number generator is common to all instances of the class SupressCrownTop
    static Uniform u;
  };

  ///\brief Reduce apicality by deleting buds in the top part of the crown
  ///\tparam TS Tree segment
  ///\tparam BUD Bud
  template <class TS,class BUD>
  class ReduceApicalityDeleteBuds{
  public:
    ///\brief Constructor
    ///\param r_height Relative height in the interval [0,1] defining the top part of the tree crown
    ///\param p Probability of a bud in the top part of the crown to die
    ReduceApicalityDeleteBuds(double r_height, double p)
      :crown_top(r_height),p_dead(p){}
    ///\brief Reduce apicality by killing buds in the upper part of the tree crown
    ///\param t Tree
    ///\note Meant to be used with std::for_each and the tree vector in GrowthLoop
    ///\sa SupressCrownTop
    void operator()(Tree<TS,BUD>* t)const;
  private:
    ///Top part of the tree crown in the interval [0,1].
    ///For example buds \f$ 0.9 \times \mathrm{tree_height}\f$ and above are in the top part of the crown.
    double crown_top;
    double p_dead;///< Probability of a bud in the top part of the crown to die
  };
  
  template<class TS, class BUD>
  TreeCompartment<TS,BUD>* SupressCrownTop<TS,BUD>::operator()(TreeCompartment<TS,BUD>* tc)const
  {
    if (BUD* b = dynamic_cast<BUD*>(tc) && (GetValue(*b,LGAstate) == ALIVE)){
      double z = GetPoint(*b).getZ();
      if (z >= crown_top*tree_height){
	double r_dead = u(1);
	if (r_dead <= p_dead){
	  SetValue(*b,LGAstate,DEAD);
	}
      }
    }
    return tc;
  }

  template<class TS, class BUD>
  SupressCrownTop<TS,BUD>& SupressCrownTop<TS,BUD>::operator=(const SupressCrownTop& sct)
  {
    tree_height = sct.tree_height;
    crown_top = sct.crown_top;
    p_dead = sct.p_dead;
    return *this;
  }

  template <class TS,class BUD>
  void ReduceApicalityDeleteBuds<TS,BUD>::operator()(Tree<TS,BUD>* t)const
  {
    Point p(0,0,0);
    FindHighestPoint<TS,BUD> fhp;
    p = Accumulate(*t,p,fhp);
    SupressCrownTop<TS,BUD> sct(p.getZ(),crown_top,p_dead);
    ForEach(*t,sct);
  }

  ///\brief Reduce apicality in the top part of the crown
  ///
  ///Reduce apicality based on tree age, segment relative height in tree crown and Gravelius order.
  ///The function of Gravelius order is implemented as ParametericCurve in a file.
  ///Tree age and segment relative height value trigger apicality reduction.
  ///\sa  LignumForest::REDUCE_APICAL_AGE LignumForest::TREE_CROWN_TOP  LignumForest::REDUCE_APICAL_FILE 
  class ReduceApicalityInScotsPineSegmentLength{
  public:
    ///\brief Constructor where apicality parameter values are from global variables
    ///\param lambda Lambda to iterate over new tree segment lengths
    ///\sa LignumForest::REDUCE_APICAL_AGE LignumForest::TREE_CROWN_TOP LignumForest::REDUCE_APICAL_FILE
    ReduceApicalityInScotsPineSegmentLength(double lambda)
      :l(lambda),tree_apical_age(LignumForest::REDUCE_APICAL_AGE),crown_top(LignumForest::TREE_CROWN_TOP),
       fapicality(LignumForest::REDUCE_APICAL_FILE),tree_height(0),tree_height_calculated(false){}
    ///\brief Expicitely give all data
    ///\param lambda Lambda to adjust segment length
    ///\param age Tree age
    ///\param top Top part of the tree crown
    ///\param apical Apicality reduction as a function of Gravelius order
    ReduceApicalityInScotsPineSegmentLength(double lambda, double age, double top, const ParametricCurve& apical)
      :tree_apical_age(age),crown_top(top),fapicality(apical),tree_height(0),tree_height_calculated(false){}
    ///\brief Segment length with decreased growth potential for reduced apicality
    ///
    ///Segment length is based on LignumForest::SetScotsPineSegmentLength Basic model
    ///with addition to reduce apicality in the tree. When and how reduced apicality is applied
    ///depends on tree age and segment position in the tree crown.
    ///\param tc Tree compartment
    ///\retval tc The tree comaprtment
    ///\sa GrowthLoop::allocationAndGrowth LignumForest::SetScotsPineSegmentLength
    TreeCompartment<ScotsPineSegment,ScotsPineBud>* operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const;
  private:
    ///Lambda to iterate over new tree segment lengths
    double l;
    ///The tree age when the relative growth reduction will start
    ///\sa LignumForest::REDUCE_APICAL_AGE
    double tree_apical_age;
    ///Top part of the crown, \p crown_top \f$ \in [0,1] \f$
    ///\sa LignumForest::TREE_CROWN_TOP
    double crown_top;
    ///Reduce relative growth potential as a function of Gravelius order, \p fapicality \f$ \in [0,1] \f$.
    ///\sa LignumForest::REDUCE_APICAL_FILE
    ParametricCurve fapicality;
    ///Heighest Z coordinate in a tree
    ///\sa Lignum::FindHighestPoint
    mutable double tree_height;
    ///Boolean value if the tree height has been calculated once \sa ReduceApicalityInScotsPineSegmentLength::operator()
    mutable bool tree_height_calculated;
  };

  TreeCompartment<ScotsPineSegment,ScotsPineBud>* ReduceApicalityInScotsPineSegmentLength::operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (Axis<ScotsPineSegment,ScotsPineBud>* axis = dynamic_cast<Axis<ScotsPineSegment,ScotsPineBud>*>(tc)){
      if (!tree_height_calculated){
	Point p(0,0,0);
	ScotsPineTree& t = dynamic_cast<ScotsPineTree&>(GetTree(*axis));
	p = Accumulate(t,p,FindHighestPoint<ScotsPineSegment,ScotsPineBud>());
	tree_height = p.getZ();
	tree_height_calculated=true;
      }
    }	  
    ///\par Steps in setting segment length
    ///The first part is as in SetScotsPineSegmentLength. The final part reduces apicality.
    else if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0.0){
	///\par Relative light
	///Mandatory relative light effect on segment length with Basic model
	///\internal
	///\snippet{lineno} ApicalDominance.h ApicalRL
	// [ApicalRL]
	double B = GetValue(GetTree(*ts),TreeQinMax);
	double qin = GetValue(*ts,LGAQin);
	double ip = qin/B;
	// [ApicalRL]
	///\endinternal
	double go =  GetValue(*ts,LGAomega); 
	double Lnew = 0.0;
	///\par Vigour index
	///For side branches the *vigour index* (Tree Physiology) is needed
	///\sa Lignum::TreePhysiologyVigourIndex
	///\internal
	///\snippet{lineno} ApicalDominance.h AVigorIndex
	// [AVigorIndex]
	double vi = GetValue(*ts,LGAvi);
	// [AVigorIndex]
	///\remark The *growth vigour* effect on segment length is a function of vigor index \f$ \mathit{f(vi)} \f$, for example
	///\f[
	///\mathit{f(vi)} = (1-a)f(\mathit{vi}) = (1-0.2)(0.15+0.85\mathit{vi}) = 0.8(0.15+0.85\mathit{vi})
	///\f]
	///\par  Growth functions
	///Functions and parameters for Basic model can be be queried from the tree:
	/// -# \p fgo: Segment length as a function of Gravelius order
	/// -# \p fip: Segment length as a function of relative light
	/// -# \p fvi: Segment length as a function of Vigor index, i.e. *growth vigor*
	/// -# \p faf: Initial foliage for Scots pine as a function of relative light
	/// -# \p flr: Segmment radius as a function of relative light
	/// -# \p getApical(): Function to set apicality [0,1] for new segments, 1 for others
	/// -# \p fapicality: Experimental addition. Apicality in older trees set for new segments in the upper part of the crown.
	///                   Set ScotsPineSegment::getApical to 1 when using \p fapicality.
	///\internal
	///\snippet{lineno} ApicalDominance.h ApicalFunc
	// [ApicalFunc]
	ParametricCurve fgo(GetFunction(GetTree(*ts),Lignum::LGMGO));
	ParametricCurve fip(GetFunction(GetTree(*ts),Lignum::LGMIP));
	const ParametricCurve& fvi = GetFunction(GetTree(*ts),LGMVI);
	const ParametricCurve& faf = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFAF);
	const ParametricCurve& flr = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFLR);
	// [ApicalFunc]
	///\endinternal
	///
	///\par Segment apicality in the Basic model
	///The value of  \p my_apical is [0,1] for new branches  and set to 1 after that.
	///Also set apicality \p reduction to 1.0 (no effect) in the tree crown apicality model for older trees.
	///\internal
	///\snippet{lineno} ApicalDominance.h ReduceApicalityBasic
	// [ReduceApicalityBasic]
	double my_apical = ts->getApical();
	double reduction = 1.0;
	// [ReduceApicalityBasic]
	///\endinternal
	//Apical reduction in upper part of the crown
	double z = GetPoint(*ts).getZ();
	if ((GetValue(GetTree(*ts),LGAage) >= tree_apical_age) && ((z / tree_height) >= crown_top)){
	  ///\par Reduce apicality
	  ///If the tree age and the segment height  in the tree crown triggers the tree crown apicality model
	  ///reduce the segment elongation defined by Gravelius order.
	  ///Also set \p my_apical in the Basic model to 1.0 (no effct)
	  ///\internal
	  ///\snippet{lineno} ApicalDominance.h LApicalReduction
	  // [LApicalReduction]
	  my_apical=1.0;
	  reduction = fapicality(go);
	  // [LApicalReduction]
	  ///\endinternal
	}
	///\par Final segment length
	///Final segment length for the Basic model with the apicality reduction in the
	///top part of the tree crown
	///\internal
	///\snippet{lineno} ApicalDominance.h LApicalBasic
	// [LApicalBasic]
	Lnew = l*reduction*my_apical*fip(ip)*fvi(vi)*fgo(go);
	// [LApicalBasic]
	///\remark Either \p reduction (apicality in the top part of the tree crown)
	///or \p my_apical (Basic model) is set to 1.0
	///\endinternal
	///
	///\par Set segment dimensions
	///Set new segment dimensions as in SetScotsPineSegmentLength
	///\internal
	///\snippet{lineno} ApicalDominance.h LApicalDimensions
	// [LApicalDimensions]
	//Set segment dimensions based on segment length
	SetValue(*ts,LGAL,Lnew);
	//Initial radius, NOTE: R = flr(ip)*Lnew, i,.e. length/radius function NOT static constant  
	SetValue(*ts,LGAR,flr(ip)*Lnew);
	//Reset previous Rh!!!!
	SetValue(*ts,LGARh,0.0);
	//Initial heartwood
	SetValue(*ts,LGARh,sqrt((GetValue(GetTree(*ts),LGPxi)*GetValue(*ts,LGAAs))/PI_VALUE));
	//SetValue(*ts,LGAWf,faf(ip)*GetValue(*ts,LGASa));
	//Needle mass now proportional to segment length
	SetValue(*ts,LGAWf,faf(ip)*Lnew);   
	//Remember the initial foliage!!
	SetValue(*ts,LGAWf0,GetValue(*ts,LGAWf));
	//Remember original sapwood area As0
	SetValue(*ts,LGAAs0,GetValue(*ts,LGAAs));
	// [LApicalDimensions]
	///\endinternal
      }//GetValue(*ts,LGAage)==0
    }//ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment>(tc)
    return tc;
  }//Function implementation end
}//namespace LignumForest
#endif
