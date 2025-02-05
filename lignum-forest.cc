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
///Use CMake to compile LignumForest. The binary after compilation is called `lignum-forest`
///With CMake all the compilation work is done in a designated compilation directory using CMakeLists.txt files.
///There are several L-system files to experiment with tree architecture. Check LSYSTEMFILE and LSYSTEMSRC variables
///in CMakeLists.txt before compilation. 
///\attention In compilation `make install` is required to copy `lignum-forest` to LignumForest working directory. 
///\deprecated The Qt `qmake` build tool is obsolete. See the README.md file for details.
///\par CMakeLists.txt
///\include CMakeLists.txt
///\page cmakefileforest

///\page runscriptforest Run LignumForest 
///The following  `run-lignum-forest.slurm` script can be used to run `lignum-forest`
///either using `bash` shell or with Slurm worload manager.
///Check the command line parameters for the `lignum-forest` binary.
///The `-metafile` option accepts glob expressions. Protect the glob expression with quotes.
///\note Remember to update the option *-hdf5* file name.
///\par run-lignum-forest.slurm
///\include run-lignum-forest.slurm
///\page runscriptforest

///\page lsystemcforest L-system experiment C
///The L-system implemented in pine-em98-branch-C.L is identical to the one in CrownDensity.
///\par Physiological age
///Introduce concept *physiological age* for buds. Each branch bud inherits the physiological
///age of the mother bud  and each time step increase the physiological age by 1. When the physiological age
///reaches  LignumForest::architecture_change_year the terminating bud genererates always exactly *two* side branches.
///Only Turn, no Roll along mother axis. These side branches have high branching angle to make denser tree crown.
///In this case branching stops if and only if tree parameters and tree functions determine so.
///\sa  `Pine::turn_branch_max` value in pine-em98-branch-C.L.
///\par pine-em98-branch-C.L
///\include pine-em98-branch-C.L
///\page lsystemcforest

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
  /// \defgroup  AMAIN The LignumForest main program
  /// @{
  /// \par The main program
  /// + Create the forest stand
  /// + HDF5 files set-up
  /// + Crowth initialization
  /// + Growth loop
  ///   + Prepare growth step
  ///   + Growth step
  ///   + Data collection
  /// + After growth tasks
  /// @}
  Sensitivity<ScotsPineSegment,ScotsPineBud> sensitivity;
  ScotsPineForest gloop;
  ran3(&LignumForest::ran3_seed);
  /// \defgroup BINIT Create the forest stand
  /// \addtogroup AMAIN
  /// @{
  /// \par Create the forest stand
  /// + Parse and resolve command line
  /// + Initialize functions
  /// + Create tree locations
  /// + Create trees
  ///
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
  // [InitForest]
  /// @}
  /// \defgroup CHDF5 HDF5 files set-up
  /// \addtogroup AMAIN
  /// @{
  /// \par HDF5 files set-up
  /// + Create HDF5 file for forest  stand data
  /// + Create HDF5 file for XML trees
  /// + Create HDF5 group for XML trees
  /// + Create HDF5 datasets for simulation configuration
  ///   + Command line
  ///   + MetaFiles
  ///   + Parameters
  ///   + Functions
  ///   + Firmament
  ///   + Initial VoxelSpace
  ///
  /// \snippet{lineno} lignum-forest.cc HDF5Init
  // [HDF5Init]
  string hdf5fname;
  ParseCommandLine(argc,argv,"-hdf5", hdf5fname);
  LGMHDF5File hdf5_trees(TREEXML_PREFIX+hdf5fname);
  hdf5_trees.createGroup(TXMLGROUP);
  LignumForest::CreateHDF5File hdf5datafile(hdf5fname,gloop.getVoxelFile(),gloop.getMetaFiles());
  hdf5datafile.createConfigurationDataSets(argc,argv);
  // [HDF5Init]
  /// @}
  /// \defgroup DGROWTHINIT Growth initialization
  /// \addtogroup AMAIN
  /// @{
  /// \par Growth initialization
  ///  + Initialize trees
  ///  + Resize HDF5 data arrays
  ///  + Initialize voxel space
  ///  + Evaluate stand variables
  ///  + Collect the initial forest data
  ///
  /// \snippet{lineno} lignum-forest.cc InitForestGrowth
  // [InitForestGrowth]
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
  // [InitForestGrowth]
  /// @}
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
    /// \defgroup EPREPARE Prepare growth step
    /// \addtogroup AMAIN
    /// @{
    /// \par Prepare growth step
    /// + Pass growth \p year to L system and to \c gloop.
    /// + Save previous year tree height for each tree.
    /// + Save the current \p year in LignumForest::GrowthLoop.
    /// + Increase \f$\xi\f$, parameter for the share of heartwood in new segments.
    /// + Growth mode change check.
    ///
    /// \sa Lignum::LGPxi
    ///
    /// \snippet{lineno} lignum-forest.cc InitGrowthLoop
    // [InitGrowthLoop]
    //Pine::L_age is for L-system 
    Pine::L_age = (double)year;     
    gloop.setHPrev();
    gloop.setYear(year);
    gloop.increaseXi(year);
    gloop.growthModeChange(year);
    // [InitGrowthLoop]
    /// @}
    ///
    /// \defgroup  FSTEP Growth step
    /// \addtogroup AMAIN
    /// @{
    /// \par Growth step
    ///  + Update foliage in voxel space and recalculate border forest
    ///  + Calculate radiation climate for trees
    ///  + Calculate photosynthesis, respiration and aging of tree compartments. 
    ///  + Create new segments
    ///  + Growth allocation
    ///  + Prune dead branches from trees
    ///  + Set radiation use efficiency in new segments (command line)
    ///    + Function of shadiness experienced by mother segment (command line argument)
    /// \sa GrowthLoop::photosynthesisRespirationTreeAging().
    ///
    /// \snippet{lineno} lignum-forest.cc NewSeg
    // [NewSeg]
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
    // RUE: radiation use efficiency
    gloop.radiationUseEfficiency();
    // [NewSeg]
    /// @}
    ///
    /// \defgroup GDATA Data collection
    /// \addtogroup AMAIN
    /// @{
    /// \par Data collection
    ///  + Evaluate stand metrics
    ///  + Collect VoxelSpace dimensions
    ///  + Collect tree data 
    ///  + Save trees in XML format in HDF5 file with write intervals
    ///
    /// \snippet{lineno} lignum-forest.cc DataCollection
    // [DataCollection] 
    gloop.evaluateStandVariables();
    // collectDataAfterGrowth collects data for HDF5 file. The 0th Year dimension
    // contains initial data. The year+1 requires because year is not yet updated
    // to denote the end of the growth
    // and the year=0 is reserved of the initial data. 
    gloop.collectDataAfterGrowth(year+1);
    //Collect VoxelSpace dimensions
    gloop.collectVoxelSpaceData(year+1,gloop.getWriteInterval());
    //Save trees as xml
    CreateTreeXMLDataSet(gloop,hdf5_trees,TXMLGROUP,gloop.getWriteInterval());
    // [DataCollection]
    /// @}
    ///
  } // End of  for(year = 0; ...)
  /// \defgroup HAFTERGROWTH After growth tasks
  /// \addtogroup AMAIN
  /// @{
  /// \par After growth tasks
  ///   + Clean up growth loop
  ///   + Collect data in HDF5 files
  ///      + Year by year, tree by tree data
  ///      + Aggregate stand data
  ///      + Aggregate center stand data
  ///      + Save HDF5 files
  ///
  /// \snippet{lineno} lignum-forest.cc AfterGrowth
  // [AfterGrowth]
  gloop.cleanUp();
  if (gloop.getNumberOfTrees() == 0){
    cout << "NO TREES AFTER GROWTH LOOP" <<endl;
  }
  cout << "GROWTH DONE " << "NUMBER OF TREES " << gloop.getNumberOfTrees() << endl;

  hdf5datafile.createDataSets(gloop);
  hdf5datafile.close();
  // [AfterGrowth]
  /// @}
  ///
  cout << "HDF5 DATA SAVED AND SIMULATION DONE" <<endl;
  return 0;
}
