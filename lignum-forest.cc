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
///\brief Declaration of a number of global variables.
///
///They are easy to add to the program but maybe should be implemented as function arguments or something else,
///see e.g. PineBudData in Pine.h \sa PineBudData
/// @{

bool is_adhoc = false;///< If function adhoc is in use: it boosts shoot growth in lower parts of the crown, see Eq. 8 in Sievanen et al. 2018 \sa adhoc
ParametricCurve adhoc("adhoc.fun");///< This function boosts shoot growth in lower parts of the crown, see Eq. 8 in Sievanen et al. 2018 \sa SetScotsPineSegmentLength
double global_hcb;///< Height of grown base for function adhoc \sa adhoc

///\defgroup spaceo Global variables to asses space occupation in shoot growth
/// @{
/// \brief Eq 10 in Sievanen et al. 2018 \sa SetScotsPineSegmentLength

Firmament dummy_firm;
/// If foliage in growth direction, Eq 10 in Sievanen et al. 2018 \sa SetScotsPineSegmentLength
VoxelSpace space_occupancy(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
			   0.1,0.1,0.1,5,5,5,dummy_firm);
///only voxelbox at the end of new Segment is checked
bool space0 = false;
///also neighboring boxes in direction of new Segment are checked
bool space1 = false;
///also all neighboring boxes are checked
bool space2 = false;
///search distance of neighboring boxes for space2
double space2_distance = 0.3;
///@}

///\defgroup spaceb Global variables to asses amount of foliage in front of a bud
/// @{
/// \brief Eq 11 in Sievanen et al. 2018

///This global function is used to asses
///foliage area density in the view cone of a bud. It affects
///the number of lateral buds created in L-system \sa pine-em98.L  \sa SetBudViewFunctor
ParametricCurve bud_view_f;
///If *bud view function* is in use. \sa bud_view_f 
bool is_bud_view_function = false;   
///@}

///These global variables have been declared in L-system and convey 
///tree age and height to L-system. \sa pine-em98.L
extern double L_age, L_H;

///Random seed for `ran3` function
int ran3_seed;

///Initial height ot trees \sa pine-em98.L
double H_0_ini;

///\defgroup  ranvar Global variables for variation
/// @{
/// \brief Variables to generate random variation between tree individuals.
/// They work in the Lindenmayer system \sa pine-em98.L

///Variation of initial heights
double  H_var_ini;
 ///Variation in the initial number of buds.
int n_buds_ini_min, n_buds_ini_max;
///Variation in the number of buds
double rel_bud;
bool bud_variation; ///< If bud variation is on \sa pine-em98.L
double branch_angle; ///< Variation in branching_angle
///@}

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
  /// \subsection varsg Growth loop variables
  /// \internal
  /// \snippet{lineno} lignum-forest.cc Vars
  // [Vars]
  Sensitivity<ScotsPineSegment,ScotsPineBud> sensitivity;                                                           
  GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > gloop;
  // [Vars]
  /// \endinternal
  
  ran3(&ran3_seed);

  /// \internal
  /// MixedForest  will  implement a  class  that  can  manage a  forest
  /// consting of several tree species.
  /// This is a quick way to test that the class compiles.
  /// \snippet{lineno} lignum-forest.cc Mf
  // [Mf] 
  // MixedForest<ScotsPineForest,ScotsPineForest> mf;
  // mf.initialize(argc,argv);
  // mf.growthLoop();
  // mf.afterGrowth();
  // [Mf]
  /// \endinternal
  
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
