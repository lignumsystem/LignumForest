#ifndef HARVEST_FOREST_STAND_H
#define HARVEST_FOREST_STAND_H

#include <mathsym.h>
#include <vector>
#include <Lignum.h>
#include <ScotsPine.h>
 
using namespace std;
using namespace cxxadt;

namespace LignumForest{
int HarvestForestStand(vector<pair<double,double> >& v,
		       vector<Tree<ScotsPineSegment,ScotsPineBud>*>& tv,
		       vector<Pine::LSystem<ScotsPineSegment,ScotsPineBud,
		       PBNAME,PineBudData>*>& lv, double r);
void RemoveDeadTrees(vector<pair<double,double> >& v,
		     vector<Tree<ScotsPineSegment,ScotsPineBud>*>& pinev,
		     vector<Pine::LSystem<ScotsPineSegment,ScotsPineBud,
		     PBNAME,PineBudData>*>& plv);
}
		     
#endif
