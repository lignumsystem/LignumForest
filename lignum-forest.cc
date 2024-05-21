/// \file lignum-forest.cc 
/// \brief The main program for the LignumForest.
///
///Main program growth loop generates now only HDF5 files for result analysis.<br>
///Growth loop steps
///+ Initialize global variables
///+ Read command line
///+ Initialize tree, L-system and global variables from command line
///+ Run the simulation
///+ Save simulation data to HDF5 files.
///
///\page cmakefileforest Compile LignumForest
///CMakeLists.txt to compile LignumForest. There are several L-system files to experiment
///with tree architecture. Check LSYSTEMFILE and LSYSTEMSRC variables in CMakeLists.txt
///before compilation. With CMake all the compilation work is done in a designated
///compilation directory. 
///\note The binary after compilation is called `lignum-forest`
///\attention Remamber to type `make install` in the compilation directory
///to copy `lignum-forest` to *LignumForest* working direcgtory.
///\deprecated The Qt `qmake` build tool is obsolete (see the README.md file for details).
///
///\include CMakeLists.txt
///\page runscriptforest Run LignumForest 
///Use the following `run-lignum-forest.sh` script to run `lignum-forest` with the latest
///CrownDensity (i.e. `crowndens`) parameters and function set. The command line for `lignum-forest` must be
///checked to be synchronised with `crowndens` before any serious simulations, most notably
///for growth mode and architectural change years.
///\note The *-metafile* option accepts regular expressions. Protect the regular expression with quotes.
///\note Remember to update the option *-hdf5* file name.
///\note In compilation `make install` is required to copy `lignum-forest` to LignumForest working directory. 
///
///\include run-lignum-forest.sh
///
///\page lsystemCforest L-system experiment C
///The L-system implemented in pine-em98-branch-C.L is identical to the one in CrownDensity. 
///Introduce concept *physiological age* for buds. Each branch bud inherits the physiological
///age of the mother bud  and each time step increase the physiological age by 1. When the physiological age
///reaches  LignumForest::architecture_change_year the terminating bud genererates *always* two side branches.
///Only Turn, no Roll along mother axis. *These side branches have high branching angle*.
///See `Pine::turn_branch_max` value in pine-em98-branch-C.L.
///In this case branching stops if and only if tree parameters and tree functions determine so.
///
///\include pine-em98-branch-C.L


//Lignum implementation
#include <Lignum.h>
#include <LignumForestGlobals.h>
#include <CreateHDF5Files.h>
#include <GrowthLoop.h> 
//Implementation of the Scotspine tree segment and bud
#include <ScotsPine.h>
#include <SomeFunctors.h>
#include <DiameterGrowth.h>
#include <Palubicki_functors.h>
#include <Space.h>


#if defined (__APPLE__) || defined(__MACOSX__)
#include <VisualFunctor.h>
//Impelements VisualizeLGMTree
#include <GLSettings.h>
#include <OpenGLUnix.h>
#include <LGMVisualization.h>
#endif
//Includes L-system, turtle graphics, GrowthLoop typedefs
#include <GLoopDefs.h>

using namespace Pine;
using namespace LignumForest;



int main(int argc, char** argv)
{ 
  /// \defgroup  lignumforest The LignumForest main program
  ///@{
  /// \par The main growth loop 
  /// + Initialize the forest stand
  /// + Set-up HDF5 file
  /// + Growth and data collection
  /// + Write data to HDF5 files
  ///@}
  ///
  Sensitivity<ScotsPineSegment,ScotsPineBud> sensitivity;                                                           
  ScotsPineForest gloop;
  ran3(&LignumForest::ran3_seed);
  ///
  /// \defgroup A_initstand Initialize forest stand
  /// \addtogroup lignumforest
  ///@{
  /// \par Initializing forest stand 
  /// + Initialize functions
  /// + Create and initialize trees to defined locations
  /// + Initialize voxel space
  /// + Resize data matrices
  /// + Collect data from initial state
  /// \internal
  /// \snippet{lineno} lignum-forest.cc InitForest
  // [InitForest]
  gloop.parseCommandLine(argc,argv);
  gloop.resolveCommandLineAttributes();
  gloop.printVariables();
  //Read and install functions
  gloop.initializeFunctions();
  //Generate tree locations either on the fly or from a file
  gloop.setTreeLocations();
  //Create trees to locations set above
  cout << "CREATE TREES" <<endl;
  gloop.createTrees();
  cout << "CREATE TREES DONE" << endl;
  gloop.printTreeLocations(0);
  //InitializeTrees reads in/sets a number of parameters and functions for each tree
  gloop.initializeTrees();
  //Resize the 3D and 2D data arrays for HDF5 file to right dimensions.
  //Simulation years and number of trees are known
  gloop.resizeTreeDataMatrix();
  gloop.initializeVoxelSpace();
  //Sets initial values of some variables in trees
  gloop.initializeGrowthLoop();
  // Evaluate stand variables before collectDataAfterGrowth
  gloop.evaluateStandVariables();
  //The 0th Year dimension is used for intial data 
  gloop.collectDataAfterGrowth(0);
  // [InitForest]
  ///@}
  /// \endinternal
  ///
  /// \defgroup B_hdf5setup Set-up HDF5 files
  /// \addtogroup lignumforest
  ///@{
  /// \par Setting up HDF5 files
  /// + Create HDF5 file for forest  stand data
  /// + Create HDF5 file for XML trees
  /// + Create HDF5 group for XML trees
  /// \internal
  /// \snippet{lineno} lignum-forest.cc HDF5INIT
  // [HDF5INIT]
  string hdf5fname;
  ParseCommandLine(argc,argv,"-hdf5", hdf5fname);
  LGMHDF5File hdf5_trees(TREEXML_PREFIX+hdf5fname);
  hdf5_trees.createGroup(TXMLGROUP);
  // [HDF5INIT]
  ///@}
  /// \endinternal
  ///
  cout << "INIT DONE" << endl;
  for(int year = 0; year < gloop.getIterations(); year++) {
    cout << "GROWTH LOOP YEAR " << year <<endl;
    if(gloop.getNumberOfTrees() < 1) {
      cout << "Number of trees left " << gloop.getNumberOfTrees() << " Stop." << endl;
      //Do not stop abruptly, continue to the end of loop and then write data
      //HDF5 output assumes there is at least one tree left 
      continue;
    }
    ///\defgroup C_datatransfer Save data from previous year
    ///\addtogroup lignumforest
    ///@{
    ///\par Save data
    ///+ Pass growth \p year to L system and to *gloop*.
    ///+ Save previous year tree height for each tree.
    ///+ Save the current \p year in LignumForest::GrowthLoop.
    ///+ Increase \f$\xi\f$, parameter for the share of heartwood in new segments.
    ///+ Growth mode change check.
    ///\sa Lignum::LGPxi.
    ///\internal
    ///\snippet{lineno} lignum-forest.cc IGLOOP
    // [IGLOOP]
    //Pine::L_age is for L-system 
    Pine::L_age = (double)year;     
    gloop.setHPrev();
    gloop.setYear(year);
    gloop.increaseXi(year);
    gloop.growthModeChange(year);
    // [IGLOOP]
    ///@}
    ///\endinternal
    ///\defgroup D_radandgrowth Radiation regime and new segments
    ///\addtogroup lignumforest
    ///@{
    ///\par The steps from the current state to a new state in the forest plot 
    ///+ Update foliage in voxel space and recalculate borderforest
    ///+ Calculate radiation climate for trees
    ///+ Calculate photosynthesis and respiration as well as aging of tree compartments. 
    ///+ Create new segments
    ///+ Growth allocation
    ///+ Prune dead branches from trees
    ///\sa GrowthLoop::photosynthesisRespirationTreeAging().
    ///\internal
    ///\snippet{lineno} lignum-forest.cc NEWSEG
    // [NEWSEG]
    gloop.setVoxelSpaceAndBorderForest();
    gloop.calculateRadiation();
    //Currently data collection in photosynthesisRespirationTreeAging
    //before new growth and tree aging
    gloop.photosynthesisRespirationTreeAging();
    gloop.createNewSegments();
    // REMOVE THESE WHEN YOU REMOVE fgomode,fipmode FROM LGMGrowthAllocator2 !!!!!!!!!!!!!!!!!!
    ParametricCurve fip_mode = GetFunction(*(gloop.getTreeVector())[0], LGMIP);
    ParametricCurve fgo_mode = GetFunction(*(gloop.getTreeVector())[0], SPFGO);

    gloop.allocationAndGrowth(fip_mode,fgo_mode);
    // Prune dead parts from the trees 
    gloop.prune();
    // [NEWSEG]
    ///@}
    ///\endinternal
    ///\defgroup E_ruegroup Radiation use efficiency
    ///\addtogroup lignumforest
    ///@{
    ///\par Set radiation use efficiency
    ///Set  radiation use efficiency in new segments as a function of shadiness
    ///experienced by mother segment.
    ///\internal
    ///\snippet{lineno} lignum-forest.cc RUE
    // [RUE]
    gloop.radiationUseEfficiency();
    // [RUE]
    ///@}
    ///\endinternal
    ///\defgroup F_datacollect Data collection
    ///\addtogroup lignumforest
    ///@{
    ///\par Data collection steps after growth
    ///+ Evaluate stand metrics 
    ///+ Collect tree data 
    ///+ Save trees in XML format in HDF5 file with write intervals 
    ///\internal
    ///\snippet{lineno} lignum-forest.cc DATACOLLECT
    // [DATACOLLECT] 
    gloop.evaluateStandVariables();
    // collectDataAfterGrowth collects data for HDF5 file. The 0th Year dimension
    // contains initial data. The year+1 requires because year is not yet updated
    // to denote the end of the growth
    // and the year=0 is reserved of the initial data. 
    gloop.collectDataAfterGrowth(year+1);
    //Save trees as xml
    CreateTreeXMLDataSet(gloop,hdf5_trees,TXMLGROUP,gloop.getWriteInterval());
    // [DATACOLLECT]
    ///@}
    ///\endinternal
  } // End of  for(year = 0; ...)
  gloop.cleanUp();
  if (gloop.getNumberOfTrees() == 0){
    cout << "NO TREES AFTER GROWTH LOOP" <<endl;
  }
  cout << "GROWTH DONE " << "NUMBER OF TREES " << gloop.getNumberOfTrees() << endl;
  ///\defgroup G_hdf5data HDF5 data creation
  ///\addtogroup lignumforest
  ///@{
  /// \par HDF5 data from growth loop data
  ///+ Year by year, tree by tree data
  ///+ Aggregate stand data
  ///+ Aggregate center stand data
  ///+ Parameters used
  ///+ Functions known in a tree
  ///+ All functions in the run directory
  ///+ Firmamment
  ///+ Initial voxel space
  ///+ Command line
  ///+ Save HDF5 files
  ///\internal
  ///\snippet{lineno} lignum-forest.cc HDF5FILECREATE
  // [HDF5FILECREATE]
  LignumForest::CreateHDF5File hdf5datafile(hdf5fname);
  hdf5datafile.createDataSets(gloop,argv,argc);
  hdf5datafile.close();
  // [HDF5FILECREATE]
  ///@}
  ///\endinternal
  cout << "HDF5 DATA SAVED AND SIMULATION DONE" <<endl;
  return 0;
}
