/// \file lignum-forest.cc 
/// \brief The main program for the LignumForest.
///
///Main program growth loop\footnote{This file gives examples of Doxygen grouping}
///generates now only HDF5 files for result analysis.<br>
///Growth loop steps
///+ Initialize global variables
///+ Read command line
///+ Initialize tree, L-system and global variables from command line
///+ Run the simulation
///+ Save simulation data to HDF5 files.
///

/// \dir ../c++adt General purpose classes and functions

/// \dir ../c++adt/include General purpose classes and functions

/// \dir ../Firmament Firmament model for radiation calculation

/// \dir ../Firmament/include Firmament model for radiation calculation

/// \dir ../stl-lignum Lignum model and generic algorithms

/// \dir ../stl-lignum/include Lignum model and generic algorithms

/// \dir ../stl-voxelspace Voxel space model

/// \dir ../stl-voxelspace/include Voxel space model and radiation calculations

/// \dir ../LEngine L-system framework

/// \dir ../LEngine/include L-system framework implementation 

/// \dir ../XMLTree  Lignum XML file format

/// \dir ../LignumForest/include Scots pine forest stand with independent trees  

/// \dir ResultAnalysis R and Python data analysis 

///\page cmakefileforest CMakeLists.txt for LignumForest
///Use CMake to compile LignumForest. The binary after compilation is called `lignum-forest`
///With CMake all the compilation work is done in a designated compilation directory using CMakeLists.txt files.
///
///There are several L-system files to experiment with tree architecture. Check LSYSTEMFILE and LSYSTEMSRC variables
///in CMakeLists.txt before compilation. 
///
///\attention In compilation `make install` is required to copy `lignum-forest` to LignumForest working directory. 
///\deprecated The Qt `qmake` build tool is obsolete. See the README.md file for details.
///\par CMakeLists.txt
///\include{lineno} CMakeLists.txt
///\page cmakefileforest

///\page runscriptforest Run LignumForest 
///The following  `run-lignum-forest.slurm` script can be used to run `lignum-forest`
///either using `bash` shell or with Slurm worload manager. Check the command line parameters
///for the `lignum-forest` binary. The `-metafile` option accepts glob expressions.
///Protect the glob expression with quotes.
///\attention Read the `ÌMPORTANT` comments before simulations
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
using namespace lignumxml;

///\defgroup AMAIN LignumForest main program
///@{
///The growth steps for the LignumForest
///\page AMAINPAGE LignumForest growth
///@}
int main(int argc, char** argv)
{
  ///\ingroup AMAIN
  ///@{
  ///\par Variables for the main growth loop
  ///+ The `ScotsPineForest gloop` is in control of the growth loop
  ///+ The `ran3` uniform random number generator used in the program
  ///\snippet{lineno} lignum-forest.cc GLoopVar
  // [GLoopVar]
  ScotsPineForest gloop;
  ran3(&LignumForest::ran3_seed);
  // [GLoopVar]
  ///\page MAINVARIABLES 1. Variables
  ///@}
  Sensitivity<ScotsPineSegment,ScotsPineBud> sensitivity;
  ///\ingroup AMAIN
  ///@{
  ///\par Steps to set-up a forest stand  
  /// + Parse and resolve command line
  /// + Initialize functions
  /// + Create tree locations
  /// + Create trees
  ///\snippet{lineno} lignum-forest.cc InitForest
  // [InitForest]
  gloop.parseCommandLine(argc,argv);
  gloop.resolveCommandLineAttributes();
  gloop.printVariables();
  //Read and install functions
  gloop.initializeFunctions();
  //Generate tree locations either on the fly or from a file
  gloop.setTreeLocations();
  //Create trees to locations set above
  cout << "Create trees" <<endl;
  gloop.createTrees();
  cout << "Create trees done" << endl;
  gloop.printTreeLocations(0);
  // [InitForest]
  ///\page CREATEFORESTSTAND 2. Forest stand
  ///@}
  ///\ingroup AMAIN
  ///@{
  ///\par Steps to set-up HDF5 files
  /// + Create HDF5 file for forest stand data
  /// + Create HDF5 datasets for simulation configuration
  ///   + Command line
  ///   + MetaFiles
  ///   + Parameters
  ///   + Functions
  ///   + Firmament
  ///   + Initial VoxelSpace
  ///   + VoxelSpace size evolution
  ///   + Create HDF5 file for XML trees
  ///   + Create HDF5 group for XML trees
  ///
  ///\snippet{lineno} lignum-forest.cc HDF5Init
  // [HDF5Init]
  string hdf5fname;
  ParseCommandLine(argc,argv,"-hdf5", hdf5fname);
  LGMHDF5File hdf5_trees(TREEXML_PREFIX+hdf5fname);
  hdf5_trees.createGroup(TXMLGROUP);
  LignumForest::CreateHDF5File hdf5datafile(hdf5fname,gloop.getVoxelFile(),gloop.getMetaFiles());
  cout << "createConfigurationDataSets" <<endl;
  hdf5datafile.createConfigurationDataSets(argc,argv);
  // [HDF5Init]
  ///\page HDF5FILES 3. HDF5 files
  ///@}
  ///\ingroup AMAIN
  ///@{
  ///\par Steps in growth initialization
  ///  + Initialize trees
  ///  + Resize HDF5 data arrays
  ///  + Initialize voxel space
  ///  + Initialize GrowthLoop::terminate_buds to check runaway branches
  ///  + Evaluate stand variables
  ///  + Collect the initial forest data
  ///  + Collect the initial voxel space dimensions
  ///\snippet{lineno} lignum-forest.cc InitForestGrowth
  // [InitForestGrowth]
  //InitializeTrees reads in/sets a number of parameters and functions for each tree
  cout << "INITIALIZE trees" << endl;
  gloop.initializeTrees();
  
  //Resize the 3D and 2D data arrays for HDF5 file to right dimensions.
  //Simulation years and number of trees are known
  gloop.resizeTreeDataMatrix();
  //VoxelSpace size ininitialization 
  gloop.initializeVoxelSpace();
  //Initialize space to check runaway branches. Create large enough space
  //inside voxel space to avoid escaped branches to expand voxel space
  gloop.initializeEscapedBuds();
  //Sets initial values of some variables in trees
  gloop.initializeGrowthLoop();
  // Evaluate stand variables before collectDataAfterGrowth
  gloop.evaluateStandVariables();
  //The 0th Year dimension is used for intial data 
  gloop.collectDataAfterGrowth(0);
  //The original voxel space
  gloop.collectVoxelSpaceData(0,gloop.getWriteInterval());
  // [InitForestGrowth]
  ///\page GROWTHINIT 4. Growth initialization
  ///@}
  /// 
  cout << "INIT DONE" << endl;
  ///\ingroup AMAIN
  ///@{
  ///\page GROWTHLOOP 5. Growth loop
  ///@}
  for(int year = 0; year < gloop.getIterations(); year++) {
    cout << "GROWTH LOOP BEGIN YEAR " << year <<endl;
    if(gloop.getNumberOfTrees() < 1) {
      cout << "Number of trees left " << gloop.getNumberOfTrees() << " Stop." << endl;
      //Do not stop abruptly, continue to the end of loop and then write data
      //HDF5 output assumes there is at least one tree left 
      continue;
    }
    /// \ingroup AMAIN
    /// @{
    /// \par Steps in preparing the  growth step
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
    ///\page PREPAREGROWTHSTEP 5.1 Prepare growth step
    /// @}
    ///
    /// \ingroup AMAIN 
    /// @{
    /// \par Steps in new growth
    ///  + Update foliage in voxel space and recalculate border forest
    ///  + Calculate radiation climate for trees
    ///  + Calculate photosynthesis, respiration and aging of tree compartments. 
    ///  + Create new segments
    ///    + Forward Qin to new segments and buds from mother segments
    ///      + Calclulate relative light (Lignum::LGAip = Lignum::LGAQin / Lignum::TreeQinMax) for terminating buds
    ///    + Set segment needle angle
    ///    + Set segment apicality
    ///    + Calculate vigour indices in tree segments and buds
    ///    + Set Lignum::TreeQinMax in the Firmament
    ///    + Calculate path length from the base of the tree to each segment
    ///    + Depending on the command line:
    ///      + Calculate space colonization model
    ///      + Calculate Extended Borchert-Honda model
    ///  + Growth allocation
    ///    + It is assumed that parameters and functions affecting segment length and diameter
    ///      can be retrieved from the segment and the Lignum tree.
    ///  + Collect data from trees that died in the allocation step
    ///  + Remove dead trees from simulation
    ///  + Check and terminate buds grown outside the Lignum::VoxelSpace
    ///  + Prune dead branches from living trees
    ///  + Depending on the command line:
    ///    + Set radiation use efficiency (RUE) in new segments
    ///    + Use function of shadiness experienced by mother segment
    /// \sa GrowthLoop::allocationAndGrowth()
    /// \sa Lignum::LGMGrowthAllocator2 and Lignum::LGMGrowthAllocator2::operator()()
    /// \sa LignumForest::SetScotsPineSegmentLength
    /// \sa LignumForest::ScotsPineDiameterGrowth2
    /// \sa  LignumForest::PartialSapwoodAreaDown
    /// \sa GrowthLoop::photosynthesisRespirationTreeAging().
    ///
    /// \snippet{lineno} lignum-forest.cc NewSeg
    // [NewSeg]
    gloop.setVoxelSpaceAndBorderForest();
    cout << "Radiation calculation" <<endl;
    gloop.calculateRadiation();
    //Currently data collection in photosynthesisRespirationTreeAging
    //before new growth and tree aging
    gloop.photosynthesisRespirationTreeAging();
    //For new segments and buds:
    //Forward Qin to new segments from mother segments
    //Forward Qin to new buds from mother segments
    //Calclulate relative light (LGAQip/TreeQinMax) for new buds
    //Set segment needle angle
    //Set segment apicality
    //Calculate vigour indices in a tree segments and buds
    //Set TreeQinMax in the Firmament
    //Calculate path length from the base of the tree to each segment
    //Depending on the command line:
    //Calculate space colonization
    //Calculate Extended Borchert-Honda model
    cout << "New growth" << endl;
    gloop.createNewSegments();
    //It assumed that parameters and functions affecting segment length and diameter
    //can be retrieved from the new segment and the Lignum tree
    gloop.allocationAndGrowth<ReduceApicalityWithFip>();
    // collectDeadTreeDataAfterGrowth collects dead tree data for HDF5 file.
    // Dead tree appears once, the year it has died and removed from simulation
    cout << "Pruning dead trees" << endl;
    gloop.collectDeadTreeDataAfterGrowth(year+1);
    // Remove dead trees from simulation
    gloop.removeDeadTreesAllOver();
    //Terminate buds grown out of VoxelSpace
    gloop.terminateEscapedBuds();
    // Prune dead parts (branches) from the trees 
    gloop.prune();
    // RUE: radiation use efficiency
    gloop.radiationUseEfficiency();
    // [NewSeg]
    ///\page GROWTHSTEP 5.2 Growth
    /// @}
    ///
    /// \ingroup AMAIN 
    /// @{
    /// \par Data collection from living trees 
    ///  + Evaluate stand metrics every year
    ///  + Collect data for each tree every year
    ///  + Collect voxelspace::VoxelSpace dimensions with write intervals
    ///  + Save trees in XML format in HDF5 file with write intervals
    ///
    /// \snippet{lineno} lignum-forest.cc DataCollection
    // [DataCollection]
    cout << "Data collection" <<endl;
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
    ///\page DATACOLLECTION 5.3 Data collection
    /// @}
    ///
    // Hard coded harvesting scheme. Define in file or command line if needed in the future
    // if (year ==  20){
    //   gloop.harvestForest(50.0);
    // }
    // if (year == 40){
    //   gloop.harvestForest(30.0);
    // }
    // if (year == 60){
    //   gloop.harvestForest(20.0);
    // }
  } // End of  for(year = 0; ...)
  /// \ingroup AMAIN
  /// @{
  /// \par Save data after growth 
  ///  + Clean up growth loop
  ///  + Create datasets for HDF5 file
  ///      + Year by year, tree by tree data
  ///      + Aggregate stand data
  ///      + Aggregate center stand data
  ///      + Voxel space dimensions data 
  ///  + Close HDF5 files for forest stand data and XML trees
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
  hdf5_trees.close();
  // [AfterGrowth]
  ///\page AFTERGROWTH 6. After growth
  /// @}
  ///
  cout << "HDF5 DATA SAVED AND SIMULATION DONE" <<endl;
  return 0;
}
