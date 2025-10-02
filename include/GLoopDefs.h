/// \file GLoopDefs.h
/// \brief GrowthLoop type alias for ScotsPineForest
#ifndef GLOOPDEFS_H
#define GLOOPDEFS_H
#include <lengine.h>
#include <ScotsPine.h>
#include <GrowthLoop.h>
#include <pinelsystem.h>
using namespace Pine;
using namespace LignumForest;

namespace LignumForest{
  ///\brief ScotsPineForest (i.e. GrowthLoop template instance) captures the growth loop of a forest stand.
  ///typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud, Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest
  typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,
  		     Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest;
}

#endif
