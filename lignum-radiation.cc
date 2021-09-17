//Include Lignum implementation
#include <Lignum.h>
#include <GrowthLoopRadiation.h> 
//Include the implementation of the tree segment and bud
#include <ScotsPine.h>

#if defined (__APPLE__) || defined(__MACOSX__)
#include <VisualFunctor.h>
//Impelements VisualizeLGMTree
#include <GLSettings.h>
#include <OpenGLUnix.h>
#include <LGMVisualization.h>
#endif
//Includes all kinds of stuff, turtle graphics etc.
#include <lengine.h>

//and for pine, see also pine9bp.L in lsys.
namespace Pine{
#include <LSystem.h>

}

bool no_compartments;
string comp_tree;

int ran3_seed = -9648383; //-23797843;

int main(int argc, char** argv)
{
  Sensitivity<ScotsPineSegment,ScotsPineBud> sensitivity;
  GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,
    Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > gloop;

  ran3(&ran3_seed);

  gloop.setVerbose(true);
  //Check and parse command line, the  command line includes switches to control
  //the simulation
  gloop.parseCommandLine(argc,argv);
  if(gloop.getOnlySelf()) {
    gloop.calculateRadiationOnlySelf();
  }
  else {
    if(gloop.getTreesFromFile()) {
      gloop.getTreesAndPositions();
    }
    else {
      gloop.setTreeLocations();
      gloop.createTrees();  //to locations set above
   
      //Write only the file about positions & trees in this case
      if(gloop.getWriteOnlyFile() && gloop.getManyTrees())
	exit(0);
    }
  }

  gloop.initializeVoxelSpace();

  gloop.setVoxelSpaceAndBorderForest();
  gloop.calculateRadiation();
}
