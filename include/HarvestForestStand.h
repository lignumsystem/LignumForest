/// \file HarvestForestStand.h
/// \brief Harvest forest stand
///
/// Harvest and remove dead trees from forest stand  
#ifndef HARVEST_FOREST_STAND_H
#define HARVEST_FOREST_STAND_H

#include <mathsym.h>
#include <vector>
#include <Lignum.h>
#include <ScotsPine.h>
 
using namespace std;
using namespace cxxadt;
namespace Pine{
  ///\relates ::LSystem
  template <class TS, class BUD, class N = LGMAD, class F=LGMdouble> class LSystem;
}

namespace LignumForest{
/// \brief Harvest forest stand
/// \param v Tree positions vector
/// \param tv Tree vector
/// \param lv Lsystem for trees vector
/// \param r Probabilty for a tree to be removed, \f$ \mathrm{r} \sim U(0,1) \f$
/// \return \f$ |v| \f$
/// \post Vectors \p v, \p tv and \p lv resized and \f$|\mathrm{v}| = |\mathrm{tv}| = |\mathrm{lv}| \f$  
int HarvestForestStand(vector<pair<double,double> >& v,
		       vector<Tree<ScotsPineSegment,ScotsPineBud>*>& tv,
		       vector<Pine::LSystem<ScotsPineSegment,ScotsPineBud,
		       PBNAME,PineBudData>*>& lv, double r);
/// \brief Remove dead trees
/// \param v Tree positions vector
/// \param pinev Tree vector
/// \param plv Lsystem for trees vector
/// \post  \f$|\mathrm{v}| = |\mathrm{pinev}| = |\mathrm{plv}| \f$  
void RemoveDeadTrees(vector<pair<double,double> >& v,
		     vector<Tree<ScotsPineSegment,ScotsPineBud>*>& pinev,
		     vector<Pine::LSystem<ScotsPineSegment,ScotsPineBud,
		     PBNAME,PineBudData>*>& plv);
}
		     
#endif
