#ifndef GLOOPDEFS_H
#define GLOOPDEFS_H
#include <lengine.h>
#include <ScotsPine.h>
#include <GrowthLoop.h>
#include <pinelsystem.h>
using namespace Pine;
using namespace LignumForest;

namespace LignumForest{
  ///typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud, Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest
  ///ScotsPineForest (i.e. GrowthLoop template instance) captures the growth loop of a forest stand.
  typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,
  		     Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest;
}

#endif
