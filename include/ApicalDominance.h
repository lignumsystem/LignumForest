#ifndef APICAL_DOMINANCE_H
#define APICAL_DOMINANCE_H
#include <Uniform.h>
/// \file ApicalDominance.h
/// \brief Experimental functors to control apical dominance in growth
///
///During later stage in development Scots pine tend to lose the strict apical
///dominance of younger trees.


namespace LignumForest{
  /// \defgroup APICALITYEXPERIMENT Apicality reduction experiment
  /// \brief Functors and functions to reduce apicality in older trees
  ///\addtogroup APICALITYEXPERIMENT
  ///@{
  ///\brief Supress growth in the top part of the crown
  ///
  ///Delete buds in top part of the tree crown. Designed to be used
  ///with Lignum::ForEach.
  ///\warning Not tested
  ///\tparam TS Tree segment
  ///\tparam BUD Bud
  template <class TS,class BUD>
  class SupressCrownTop{
  public:
    ///\brief Constructor
    ///\param t_height Tree height defined as the highest point in the tree
    ///\param r_height Relative height in the interval [0,1] defining the top part of the crown
    ///\param p Probability of a bud in the top part of the crown to die
    ///\pre \f$ 0 \leq \f$ \p r_height \f$ \leq 1 \f$
    ///\post The sequence of pseudo random numbers \p u initialized 
    SupressCrownTop(double t_height,double r_height,double p)
      :tree_height(t_height),crown_apical_top(r_height),p_dead(p){}
    ///\brief Copy constructor
    ///\param sct SupressCrownTop instance
    ///\pre SupressCrownTop::u in \p sct has been initialized.
    SupressCrownTop(const SupressCrownTop& sct)
      :tree_height(sct.tree_height),crown_apical_top(sct.crown_apical_top),p_dead(sct.p_dead){}
    ///\brief Assignment
    ///\param sct SupressCrownTop instance
    ///\return self  
    SupressCrownTop& operator=(const SupressCrownTop& sct);
    ///\brief Supress tree growth
    ///
    ///Set a bud DEAD with a probability SupressCrownTop::p_dead in the top part of the crown
    ///\param tc Tree compartment
    ///\retval tc Tree compartment 
    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;
  private:
    double tree_height;///< Highest z-coordinate in the tree
    ///Top part in the interval [0,1] of the tree crown.
    ///For example buds \f$ 0.9 \times \mathrm{tree_height}\f$ and above are in the top part of the crown.
    double crown_apical_top;
    double p_dead;///< Probability of a bud in the top part of the crown to die
    ///Uniformly distributed pseudo random number generator
    ///\note As being static the  random number generator is common to all instances of the class SupressCrownTop
    static Uniform u;
  };
  ///\brief Static data member definition and random number sequence initialization 
  template <class TS, class BUD>
  Uniform SupressCrownTop<TS,BUD>::u(-1);
  
  ///\brief Reduce apicality by deleting buds in the top part of the crown.
  ///
  ///Applied for the trees in tree vector in GrowthLoop with std::for_each.
  ///\warning Not tested
  ///\tparam TS Tree segment
  ///\tparam BUD Bud
  template <class TS,class BUD>
  class ReduceApicalityDeleteBuds{
  public:
    ///\brief Constructor
    ///\param r_height Relative height in the interval [0,1] defining the top part of the tree crown
    ///\param p Probability of a bud in the top part of the crown to die
    ReduceApicalityDeleteBuds(double r_height, double p)
      :crown_apical_top(r_height),p_dead(p){}
    ///\brief Reduce apicality by killing buds in the upper part of the tree crown.
    ///
    ///Find the highest point in the tree and use SupressCrownTop to kill buds.
    ///\note Designed to be used with std::for_each and the tree vector in GrowthLoop
    ///\param t Tree
    ///\sa SupressCrownTop
    void operator()(Tree<TS,BUD>* t)const;
  private:
    ///Top part of the tree crown in the interval [0,1].
    ///For example buds \f$ 0.9 \times \mathrm{tree_height}\f$ and above are in the top part of the crown.
    double crown_apical_top;
    double p_dead;///< Probability of a bud in the top part of the crown to die
  };

  ///\brief Reduce apicality based on Gravelius order.
  ///
  ///Reduce apicality in the top part of the tree crown based on the tree age,
  ///segment relative height in the tree crown and segment Gravelius order.
  ///The function of Gravelius order should be  implemented as ParametericCurve in a file.
  ///Tree age and segment relative height in the tree crown trigger apicality reduction.
  ///\warning Not tested
  ///\sa  LignumForest::REDUCE_APICAL_AGE LignumForest::TREE_CROWN_TOP  LignumForest::REDUCE_APICAL_FILE
  ///\sa ReduceApicalityWithFip
  ///\todo Not tested 
  class ReduceApicalityWithGraveliusOrder{
  public:
    ///\brief Constructor where apicality parameter values are set from global variables. 
    ///\param lambda Lambda to iterate over new tree segment lengths
    ///\attention This is mandatory constructor signature for GrowthLoop::allocationAndGrowth
    ///template parameter
    ///\sa LignumForest::REDUCE_APICAL_AGE LignumForest::TREE_CROWN_TOP LignumForest::REDUCE_APICAL_FILE
    ReduceApicalityWithGraveliusOrder(double lambda)
      :l(lambda),tree_apical_age(LignumForest::REDUCE_APICAL_AGE),crown_apical_top(LignumForest::TREE_CROWN_TOP),
       fapicality(LignumForest::REDUCE_APICAL_FILE),tree_height(0),tree_height_calculated(false){}
    ///\brief Expicitely give all data
    ///\param lambda Lambda to adjust segment length
    ///\param age Tree age
    ///\param top Top part of the tree crown
    ///\param apical Apicality reduction as a function of Gravelius order
    ReduceApicalityWithGraveliusOrder(double lambda, double age, double top, const ParametricCurve& apical)
      :tree_apical_age(age),crown_apical_top(top),fapicality(apical),tree_height(0),tree_height_calculated(false){}
    ///\brief Segment length with decreased growth potential for reduced apicality
    ///
    ///Segment length is based on LignumForest::SetScotsPineSegmentLength Basic model
    ///with apicality reduction according the segment Gravelius order.
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
    ///Top part of the crown, \p crown_apical_top \f$ \in [0,1] \f$
    ///\sa LignumForest::TREE_CROWN_TOP
    double crown_apical_top;
    ///Reduce relative growth potential as a function of Gravelius order, \p fapicality \f$ \in [0,1] \f$.
    ///\sa LignumForest::REDUCE_APICAL_FILE
    ParametricCurve fapicality;
    ///Heighest Z coordinate in a tree
    ///\sa Lignum::FindHighestPoint
    mutable double tree_height;
    ///Boolean value if the tree height has been calculated once \sa ReduceApicalityWithGraveliusOrder::operator()
    mutable bool tree_height_calculated;
  };
  
  ///\brief Reduce apicality in a tree based on relative light position
  ///
  ///Reduce apicality based on tree age, segment relative height in tree crown and segment relative light position.
  ///The function of relative light is implemented as ParametericCurve in a file.
  ///Tree age and segment relative height value trigger apicality reduction.
  ///\sa  LignumForest::REDUCE_APICAL_AGE LignumForest::TREE_CROWN_TOP  LignumForest::REDUCE_APICAL_FILE
  ///\sa ReduceApicalityWithGraveliusOrder
  class ReduceApicalityWithFip{
  public:
    ///\brief Constructor, apicality parameter values are set from global variables
    ///\param lambda Lambda to iterate over new tree segment lengths
    ///\attention This is mandatory constructor signature for GrowthLoop::allocationAndGrowth
    ///template parameter
    ///\sa LignumForest::REDUCE_APICAL_AGE LignumForest::TREE_CROWN_TOP LignumForest::REDUCE_APICAL_FILE
    ReduceApicalityWithFip(double lambda)
      :l(lambda),tree_apical_age(LignumForest::REDUCE_APICAL_AGE),crown_apical_top(LignumForest::TREE_CROWN_TOP),
       fapicality(LignumForest::REDUCE_APICAL_FILE){}
    ///\brief Constructor, expicitely give all data
    ///\param lambda Lambda to adjust segment length
    ///\param age Tree age
    ///\param top Top part of the tree crown
    ///\param apicality Apicality reduction as a function of Gravelius order
    ReduceApicalityWithFip(double lambda, double age, double top, const ParametricCurve& apicality)
      :tree_apical_age(age),crown_apical_top(top),fapicality(apicality){}
    ///\brief Copy constructor
    ///\param ra Existing object 
    ReduceApicalityWithFip(const ReduceApicalityWithFip& ra)
      :l(ra.l),tree_apical_age(ra.tree_apical_age),crown_apical_top(ra.crown_apical_top),
       fapicality(ra.fapicality){}
    ///\brief Assignment
    ///\param ra Existing object
    ReduceApicalityWithFip& operator=(const ReduceApicalityWithFip& ra){
      l=ra.l; tree_apical_age=ra.tree_apical_age; crown_apical_top=ra.crown_apical_top;
      fapicality=ra.fapicality; tree_height=ra.tree_height; tree_id = ra.tree_id;
      return *this;
    }
    ///\brief Segment length with decreased growth potential for reduced apicality
    ///
    ///Segment length is based on LignumForest::SetScotsPineSegmentLength Basic model
    ///with addition to reduce apicality in the tree. When and how reduced apicality is applied
    ///depends on tree age and segment position in the tree crown.
    ///\pre Tree height is calculated once.
    ///\param tc Tree compartment
    ///\retval tc The tree comaprtment
    ///\sa GrowthLoop::allocationAndGrowth()
    TreeCompartment<ScotsPineSegment,ScotsPineBud>* operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const;
  private:
    ///\brief Lambda to iterate over new tree segment lengths in ReduceApicalityWithFip::operator() 
    double l;
    ///\brief The tree age when the relative growth reduction will start
    ///\sa LignumForest::REDUCE_APICAL_AGE
    double tree_apical_age;
    ///\brief Top part of the crown, \p crown_apical_top \f$ \in [0,1] \f$
    ///\sa LignumForest::TREE_CROWN_TOP
    double crown_apical_top;
    ///\brief Reduce relative growth potential as a function of Gravelius order, \p fapicality \f$ \in [0,1] \f$.
    ///\sa LignumForest::REDUCE_APICAL_FILE
    ParametricCurve fapicality;
    ///\brief Heighest Z coordinate among living living buds in a tree
    ///\sa Lignum::FindHighestPoint
    static double tree_height;
    ///\brief Lignum::TreeId identification tag
    static double tree_id;
  };
  ///\brief Definition: the initial previous \p tree_height is nset to 0.0
  double ReduceApicalityWithFip::tree_height = 0.0;
  ///\brief Definition: the initial previous \p tree_id is set to -1.0
  double ReduceApicalityWithFip::tree_id = -1.0;

  ///\brief Kill buds after growth allocation
  ///
  ///Kill buds not able to create new segments after growth allocation.
  ///Kill buds in the top part of the tree to reduce apicality. For the latter check:
  /// + Tree age
  /// + Bud relative height in the tree
  /// + Bud has Gravelius order 1 or 2
  ///
  ///The functor for Lignum::PropagateUp() in GrowthLoop::allocationAndGrowth() and passed as template parameter.
  ///
  ///\sa LignumForest::REDUCE_APICAL_AGE
  ///\sa LignumForest::H_REL_APICAL_BUD
  ///\sa LignumForest::P_DEAD_APICAL_BUD_GO1
  ///\sa LignumForest::P_DEAD_APICAL_BUD_GO2
  ///\sa GrowthLoop::allocationAndGrowth()
  class KillBudsReduceApicality{
  public:
    ///\brief Constructor initialized with global variables to control apicality
    KillBudsReduceApicality()
      :tree_id(-1.0),tree_height(0.0),tree_apical_age(LignumForest::REDUCE_APICAL_AGE),bud_apical_rel_h(LignumForest::H_REL_APICAL_BUD),
       p_bud_dead_go1(LignumForest::P_DEAD_APICAL_BUD_GO1),p_bud_dead_go2(LignumForest::P_DEAD_APICAL_BUD_GO2){}
    ///\brief Copy constructor
    ///\param  kbra KillBudsReduceApicality object
    KillBudsReduceApicality(const KillBudsReduceApicality& kbra)
      :tree_id(kbra.tree_id),tree_height(kbra.tree_height),tree_apical_age(kbra.tree_apical_age),
       bud_apical_rel_h(kbra.bud_apical_rel_h),p_bud_dead_go1(kbra.p_bud_dead_go1),
       p_bud_dead_go2(kbra.p_bud_dead_go2){}
    ///\brief Assignment
    ///\param kbra KillBudsReduceApicality object
    ///\retval *this The object itself
    KillBudsReduceApicality& operator=(const KillBudsReduceApicality& kbra){
      tree_id = kbra.tree_id;
      tree_height = kbra.tree_height;
      tree_apical_age = kbra.tree_apical_age;
      bud_apical_rel_h = kbra.bud_apical_rel_h;
      p_bud_dead_go1 = kbra.p_bud_dead_go1;
      p_bud_dead_go2 = kbra.p_bud_dead_go2;
      return *this;
    }
    ///\brief Kill buds after growth
    ///
    ///Kill buds not able to create new segments and kill buds in the top part of the tree crown
    ///to reduce apicality. For the latter tree age, bud relative height and bud gravelius order
    ///are checked. 
    ///\param kill_bud Flag to denote if a segment has created a segment (\p kill = \e false) or not (\p kill = \e true)
    ///\param tc Tree compartment
    ///\retval kill_bud Value changed by Tree segment 
    bool operator()(bool& kill_bud, TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const;
  private:
    mutable double tree_id; ///< Lignum::TreeId identification tag
    mutable double tree_height; ///< Tree height
    double tree_apical_age;///< Tree age when the apical reduction is applied
    double bud_apical_rel_h;///< Bud relative height in tree crown when the apical reduction is applied
    double p_bud_dead_go1;///< Bud Gravelius order 1 probability to die to reduce apicality
    double p_bud_dead_go2;///< Bud Gravelius order 2 probability to die to reduce apicality
    static Uniform u;///< Random number generator for uniform distribution
  };
  ///\brief Initialize the random number generator 
  Uniform KillBudsReduceApicality::u(-1.0);
  
  template<class TS, class BUD>
  TreeCompartment<TS,BUD>* SupressCrownTop<TS,BUD>::operator()(TreeCompartment<TS,BUD>* tc)const
  {
    if (BUD* b = dynamic_cast<BUD*>(tc)){
      if (GetValue(*b,LGAstate) == ALIVE){
	double z = GetPoint(*b).getZ();
	if (z >= crown_apical_top*tree_height){
	  double r_dead = u(1);
	  if (r_dead <= p_dead){
	    SetValue(*b,LGAstate,DEAD);
	  }
	}
      }
    }
    return tc;
  }

  template<class TS, class BUD>
  SupressCrownTop<TS,BUD>& SupressCrownTop<TS,BUD>::operator=(const SupressCrownTop& sct)
  {
    tree_height = sct.tree_height;
    crown_apical_top = sct.crown_apical_top;
    p_dead = sct.p_dead;
    return *this;
  }

  template <class TS,class BUD>
  void ReduceApicalityDeleteBuds<TS,BUD>::operator()(Tree<TS,BUD>* t)const
  {
    Point p(0,0,0);
    FindHighestPoint<TS,BUD> fhp;
    p = Accumulate(*t,p,fhp);
    SupressCrownTop<TS,BUD> sct(p.getZ(),crown_apical_top,p_dead);
    ForEach(*t,sct);
  }


  TreeCompartment<ScotsPineSegment,ScotsPineBud>* ReduceApicalityWithGraveliusOrder::operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
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
	/// -# \p getApical(): Function to set apicality [0,1] for new segments, 1 for others.
	/// -# \p fapicality: Experimental addition. Apicality in older trees set for new segments in the upper part of the crown.
	///                   
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
	///The value of  \p my_apical is in [0,1] for new branches  and set to 1 after that.
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
	if ((GetValue(GetTree(*ts),LGAage) >= tree_apical_age) && (z >= crown_apical_top*tree_height)){
	  ///\par Reduce apicality
	  ///If the tree age and the segment height  in the tree crown triggers the tree crown apicality model
	  ///reduce the segment elongation defined by \p fapicality(go), where \p go is Gravelius order.
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
	//Initial radius, NOTE: R = flr(ip)*Lnew, i.e. length/radius is a function, not static constant parameter value  
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

  TreeCompartment<ScotsPineSegment,ScotsPineBud>* ReduceApicalityWithFip::operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    ///\par Calculate tree height once
    ///Save tree height for the growth iteration and set tree id to current tree id. Tree height
    ///is calculated once per growth allocation. 
    if (Axis<ScotsPineSegment,ScotsPineBud>* axis = dynamic_cast<Axis<ScotsPineSegment,ScotsPineBud>*>(tc)){
      if ((GetValue(GetTree(*axis),TreeId)) != tree_id){
	///\internal
	///\snippet{lineno} ApicalDominance.h ApicalTH
	// [ApicalTH]
	Point p(0,0,0);
	ScotsPineTree& t = dynamic_cast<ScotsPineTree&>(GetTree(*axis));
	p = Accumulate(t,p,FindHighestPoint<ScotsPineSegment,ScotsPineBud>());
	tree_height = p.getZ();
	tree_id = GetValue(GetTree(*axis),TreeId);
	// [ApicalTH]
	///\endinternal
      }
    } 
    ///\par Steps in setting segment length
    ///The first part is as in the Basic model in SetScotsPineSegmentLength.
    ///The additinional part reduces apicality when tree hase reach \p tree_apical_age
    else if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0.0){
	///\par Relative light
	///Mandatory relative light effect on segment length with Basic model
	///\internal
	///\snippet{lineno} ApicalDominance.h ApicalRLFip
	// [ApicalRLFip]
	double B = GetValue(GetTree(*ts),TreeQinMax);
	double qin = GetValue(*ts,LGAQin);
	double ip = qin/B;
	// [ApicalRLFip]
	///\endinternal
	double go =  GetValue(*ts,LGAomega); 
	double Lnew = 0.0;
	///\par Vigour index
	///For side branches the *vigour index* (Tree Physiology) is needed
	///\sa Lignum::TreePhysiologyVigourIndex
	///\internal
	///\snippet{lineno} ApicalDominance.h AVigorIndexFip
	// [AVigorIndexFip]
	double vi = GetValue(*ts,LGAvi);
	// [AVigorIndexFip]
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
	///\snippet{lineno} ApicalDominance.h ApicalFuncFip
	// [ApicalFuncFip]
	ParametricCurve fgo(GetFunction(GetTree(*ts),Lignum::LGMGO));
	ParametricCurve fip(GetFunction(GetTree(*ts),Lignum::LGMIP));
	const ParametricCurve& fvi = GetFunction(GetTree(*ts),LGMVI);
	const ParametricCurve& faf = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFAF);
	const ParametricCurve& flr = GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPFLR);
	// [ApicalFuncFip]
	///\endinternal
	///
	///\par Segment apicality in the Basic model
	///The value of  \p my_apical is [0,1] for new branches  and set to 1 after that.
	///Also set apicality \p reduction to 1.0 (no effect) in the tree crown apicality model for older trees.
	///\internal
	///\snippet{lineno} ApicalDominance.h ReduceApicalityBasicFip
	// [ReduceApicalityBasicFip]
	double my_apical = ts->getApical();
	double reduction = 1.0;
	// [ReduceApicalityBasicFip]
	///\endinternal
	//Apical reduction in upper part of the crown
	double z = GetPoint(*ts).getZ();
	if ((GetValue(GetTree(*ts),LGAage) >= tree_apical_age) && (z >= crown_apical_top*tree_height)){
	  ///\par Reduce apicality
	  ///If the tree age and the segment height  in the tree crown triggers the tree crown apicality model
	  ///reduce the segment elongation defined by relative light position. Note that \p my_apical is reset to 1.0.
	  ///\internal
	  ///\snippet{lineno} ApicalDominance.h LApicalReductionFip
	  // [LApicalReductionFip]
	  reduction = fapicality(ip);
	  my_apical = 1.0;
	  // [LApicalReductionFip]
	  ///\endinternal
	}
	///\par Final segment length
	///Final segment length for the Basic model with the apicality reduction in the
	///top part of the tree crown
	///\internal
	///\snippet{lineno} ApicalDominance.h LApicalBasicFip
	// [LApicalBasicFip]
	Lnew = l*reduction*my_apical*fip(ip)*fvi(vi)*fgo(go);
	// [LApicalBasicFip]
	///\endinternal
	///
	///\par Set segment dimensions
	///Set new segment dimensions as in SetScotsPineSegmentLength
	///\internal
	///\snippet{lineno} ApicalDominance.h LApicalDimensionsFip
	// [LApicalDimensionsFip]
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
	// [LApicalDimensionsFip]
	///\endinternal
      }//GetValue(*ts,LGAage)==0
    }//ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment>(tc)
    return tc;
  }//Function implementation end

  bool KillBudsReduceApicality::operator()(bool& kill_bud, TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  { ///\par Calculate highest point
    /// + Axis: calculate the highest point in tree nonce. Tree id tag allows the check.
    /// + Tree segment: set the \p kill_bud flag \e true if the segment could not be created.
    /// + Bud: Kill buds denoted by \p kill_bud. Kill buds meeting criteria for apical control.
    ///\internal
    ///\snippet{lineno} ApicalDominance.h ApicalControlDeleteBuds
    // [ApicalControlDeleteBuds]
    if (Axis<ScotsPineSegment,ScotsPineBud>* axis = dynamic_cast<Axis<ScotsPineSegment,ScotsPineBud>*>(tc)){
      //Calculate the highest point
      if ((GetValue(GetTree(*axis),TreeId)) != tree_id){
	Point p(0,0,0);
	ScotsPineTree& t = dynamic_cast<ScotsPineTree&>(GetTree(*axis));
	p = Accumulate(t,p,FindHighestPoint<ScotsPineSegment,ScotsPineBud>());
	tree_height = p.getZ();
	tree_id = GetValue(GetTree(*axis),TreeId);
      }
    }
    else if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      //Set the value for kill_bud tag
      METER l = GetValue(*ts,LGAL);
      if (l <= R_EPSILON){
	kill_bud = true;
      }
    }
    else if (ScotsPineBud* b = dynamic_cast<ScotsPineBud*>(tc)){
      //Kill buds 
      if (kill_bud == true){
        SetValue(*b,LGAstate,DEAD);
      }
      //Apicality control
      else{
	double age = GetValue(GetTree(*b),LGAage);
	Point p = GetPoint(*b);
	double z = p.getZ();
	if ((age >= tree_apical_age) && (z >= bud_apical_rel_h*tree_height)){
	  double go = GetValue(*b,LGAomega);
	  double p_dead = u(1);
	  if ((go == 1) && (p_dead >= 0.0) && (p_dead < p_bud_dead_go1)){
	    SetValue(*b,LGAstate,DEAD);
	  }
	  if ((go == 2) && (p_dead >= 0.0) && (p_dead < p_bud_dead_go2)){
	    SetValue(*b,LGAstate,DEAD);
	  }
	}
      }
    }
    // [ApicalControlDeleteBuds]
    ///\endinternal
    return kill_bud;
  }
 ///@} End group APICALITYEXPERIMENT
}//namespace LignumForest
#endif
