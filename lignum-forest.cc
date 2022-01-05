/// \file 
/// \brief The main program for the LignumForest.

///Lignum implementation
#include <Lignum.h>
#include <GrowthLoop.h> 
///Implementation of the Scotspine tree segment and bud
#include <ScotsPine.h>

//#include <MixedForest.h>

#if defined (__APPLE__) || defined(__MACOSX__)
#include <VisualFunctor.h>
//Impelements VisualizeLGMTree
#include <GLSettings.h>
#include <OpenGLUnix.h>
#include <LGMVisualization.h>
#endif
///Includes L-system, turtle graphics etc.
#include <lengine.h>

///L-system for pine and for pine, see also pine9bp.L in lsys.
namespace Pine{
#include <LSystem.h>
}

#include <SomeFunctors.h>
#include <DiameterGrowth.h>
//#include <RadiationCrownDens.h>
#include <Palubicki_functors.h>
#include <Space.h>

///\defgroup mainglobals Global variables
/// @{
///Declaration of a number of global variables -- they are easy to add to
///the program but maybe should be made function arguments or ...

///If growth is promoted in lower parts of crown
bool is_adhoc = false;
///Increase as a function rel. dist from crown b.
ParametricCurve adhoc("adhoc.fun");
///Height of grown base to SetScotsPineSegmentLength (if is_adhoc)
double global_hcb;

///Space colonialization options for SetScotsPineSegmentLength (in ScotsPine.h)
///roughly: a shoot can grow only if there is free space around it
///only voxelbox at the end of new Segment is checked
bool space0 = false;
///also neighboring boxes in dierction of new Segment are checked
bool space1 = false;
///also all neighboring boxes are checked
bool space2 = false;
///search distance of neighboring boxes for space2
double space2_distance = 0.3;

///Dummy firmament for space_occupancy
Firmament dummy_firm;
///VoxelSpace for space occupancy in case space colonialization
VoxelSpace space_occupancy(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
			   0.1,0.1,0.1,5,5,5,dummy_firm);

///This global variable conveys the Bud View Function to L-system
ParametricCurve bud_view_f;
bool is_bud_view_function = false;   ///if it is in use 


///These global variables have been declared in pine-em98.L and convey 
///tree age and height to L-system
extern double L_age, L_H;

/// Random seed is initialized in GrowthloopI.h
int ran3_seed;

///These variables are used to generate random variation between tree individuals
/// in the Lindenmayer systen (*.L file, as externals).
///For variation of initial heights
double H_0_ini, H_var_ini;
/// and number of buds (in .L file)
int n_buds_ini_min, n_buds_ini_max;
///Variation in no. buds
double rel_bud;
///If bud variation is on
bool bud_variation;
///For variation of branching_angle
double branch_angle;
///@}

/// \typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud, Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest
/// ScotsPineForest (i.e. GrowthLoop template instance) captures the growth loop of a forest stand.
typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,
		   Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest;

/// \defgroup groupmain Main program for growth loop
/// @{
/// Content of the main program to run Growth loop for LignumForest.
/// @}

/// \ingroup groupmain
/// \fn int main(int argc, char** argv)
/// \brief Main function for growth.
/// Check and parse command line. The  command line includes switches to control
/// the simulation.
/// \param argc An integer counter number of command line arguments.
/// \param argv Vector of command line argument strings.
int main(int argc, char** argv)
{

  /// \subsection variables Variables for the growth loop
  // \var Sensitivity sensitivity
  /// \snippet{lineno} lignum-forest.cc Vars
  /// \internal
  // [Vars]
  Sensitivity<ScotsPineSegment,ScotsPineBud> sensitivity;                                                           
  GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > gloop;
  // [Vars]
  /// \endinternal
  
  ran3(&ran3_seed);

  //MixedForest  will  implement a  class  that  can  manage a  forest
  //consting of several tree species
  //This is a quick way to test that the class compiles 
  //MixedForest<ScotsPineForest,ScotsPineForest> mf;
  //mf.initialize(argc,argv);
  //mf.growthLoop();
  //mf.afterGrowth();
  
  /// \subsection initf  Initialize forest
  /// \snippet{lineno} lignum-forest.cc InitForest
  /// \internal
  // [InitForest]
  gloop.parseCommandLine(argc,argv);
  gloop.resolveCommandLineAttributes();
  gloop.printVariables();
  gloop.initializeFunctions();
  gloop.setTreeLocations();
  gloop.createTrees();  //to locations set above
  gloop.printTreeLocations(0);
  gloop.initializeTrees();
  gloop.initializeVoxelSpace();
  gloop.initializeGrowthLoop();
  // [InitForest]
  /// \endinternal
  
  // gloop.growthLoop();

  /// \subsection growthloop Growth loop
  /// \snippet{lineno} lignum-forest.cc GLoop
  /// \internal
  // [GLoop]
  int year;
  for(year = 0; year < gloop.getIterations(); year++) {
    if(gloop.getNumberOfTrees() < 1) {
      cout << "No trees left. Stop." << endl;
      exit(0);
    }
    L_age = (double)year;     //This is for L-system and dangerous
    gloop.setHPrev();
    gloop.setYear(year);
    gloop.evaluateStandVariables();
    gloop.setVoxelSpaceAndBorderForest();
    gloop.calculateRadiation();
    gloop.increaseXi(year);
    gloop.photosynthesisAndRespiration();
    gloop.createNewSegments();
    gloop.allocationAndGrowth();
    gloop.output();   //does this if gloop.write_output == true
    //Prune dead parts from the trees 
    gloop.prune();
  } // End of  for(year = 0; ...
  // [GLoop]
  /// \endinternal
  
  /// \subsection aftegrowth After growth
  /// \snippet{lineno} lignum-forest.cc AGrowth
  /// \internal
  // [AGrowth]
  gloop.cleanUp();
  gloop.printSegmentQin();
  gloop.printBranchMeans();
  gloop.printTreeLocations(gloop.getIterations());
  gloop.printVoxelObjectLocations("VoxelObjectLocations.txt");
  gloop.writeTreeToXMLFile(gloop.getTargetTree(),GetValue(gloop.getTargetTree(),LGAage),1);
  gloop.writeFip(gloop.getTargetTree(),1);
  gloop.writeBranchInformation(gloop.getTargetTree(),"BranchInformation.dat");
  gloop.writeProductionBalance(gloop.getTargetTree(),"ProductionBalance.dat");
  // [AGrowth]
  /// \endinternal
  
#if defined (__APPLE__) || defined(__MACOSX__)
  if(CheckCommandLine(argc,argv,"-viz")) {
    LGMVisualization viz;
    viz.InitVisualization(argc,argv);
    // textures 512x512
    viz.AddCfTree(gloop.getTargetTree(), "Manty.bmp", "neulaset5.tga");
    float th = (float)GetValue(gloop.getTargetTree(),LGAH);
    cout << th << endl;
    //viz.ResetCameraPosition(th);
    viz.SetMode(SOLID);
    //viz.ResetCameraPosition(GetValue(gloop.getTargetTree(),LGAH));
    viz.StartVisualization();
  }
#endif
  return 0;
}
