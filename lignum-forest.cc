/// \file lignum-forest.cc 
/// \brief The main program for the LignumForest.

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
///before compilation.
///\note The binary is called `lignum-forest`
///\note The Qt `qmake` build tool is obsolete (see the README file for details).
///
///\include CMakeLists.txt
///\page runscriptforest Run LignumForest 
///Use the following `run-lignum-forest.sh` script to run `lignum-forest` with the latest
///CrownDensity (i.e. `crowndens`) parameters and function set. The command line for `lignum-forest` must be
///checked to be synchronised with `crowndens` before any serious simulations, most notably
///for growth mode and architectural change years. These two are important results (lessons learned)
///from successful `crowndens` simulations.
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
///Includes L-system, turtle graphics etc.
#include <lengine.h>

///L-system for pine and for pine, see also pine9bp.L in lsys.
namespace Pine{
#include <LSystem.h>
}



using namespace Pine;
using namespace LignumForest;
using namespace LignumForest;

namespace LignumForest{
  /// \typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud, Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest
  /// ScotsPineForest (i.e. GrowthLoop template instance) captures the growth loop of a forest stand.
  typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,
		     Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest;
}
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
  /// \par Growth loop variables
  /// \snippet{lineno} lignum-forest.cc Vars
  /// \internal
  // [Vars]
  Sensitivity<ScotsPineSegment,ScotsPineBud> sensitivity;                                                           
  ScotsPineForest gloop;
  // [Vars]
  /// \endinternal
  
  ran3(&LignumForest::ran3_seed);
  
  /// \par Initialize forest
  /// \snippet{lineno} lignum-forest.cc InitForest
  /// \internal
  // [InitForest]
  gloop.parseCommandLine(argc,argv);
  gloop.resolveCommandLineAttributes();
  gloop.printVariables();
  gloop.initializeFunctions();   //Reads in some functions from files
  gloop.setTreeLocations();
  gloop.createTrees();  //to locations set above
  gloop.printTreeLocations(0);
  //InitializeTrees reads in/sets a number of parameters and functions for each tree
  gloop.initializeTrees();
  //Resize the 3D and 2D data arrays for HDF5 file to right dimensions
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
  /// \endinternal
  /// \par Intialize  HDF5 content
  /// \snippet{lineno} lignum-forest.cc HDF5init
  /// \internal
  // [HDF5init]
  string hdf5fname;
  ParseCommandLine(argc,argv,"-hdf5", hdf5fname);
  LGMHDF5File hdf5_file(hdf5fname);
  LGMHDF5File hdf5_trees(TREEXML_PREFIX+hdf5fname);
  hdf5_file.createGroup(PGROUP);
  hdf5_file.createGroup(TFGROUP);
  hdf5_file.createGroup(AFGROUP);
  hdf5_trees.createGroup(TXMLGROUP);
  // [HDF5init]
  /// \endinternal
  /// \par Growth loop
  /// \snippet{lineno} lignum-forest.cc GLoop
  /// \internal
  // [GLoop]
  int year;
  for(year = 0; year < gloop.getIterations(); year++) {
    cout << "GROWTH LOOP YEAR " << year <<endl;
    if(gloop.getNumberOfTrees() < 5) {
      cout << "Number of trees left " << gloop.getNumberOfTrees() << " Stop." << endl;
      //Do not stop abruptly, continue to the end of loop and then write data
      //HDF5 output assumes there is at least one tree left 
      continue;
    }
    L_age = (double)year;     //This is for L-system and dangerous
    gloop.setHPrev();
    gloop.setYear(year);
    
    gloop.setVoxelSpaceAndBorderForest();
    gloop.calculateRadiation();
    gloop.increaseXi(year);
    // Currently data collection before new growth and tree aging
    gloop.photosynthesisAndRespiration();
    gloop.createNewSegments();

    // REMOVE THESE WHEN YOU REMOVE fgomode,fipmode FROM LGMGrowthAllocator2 !!!!!!!!!!!!!!!!!!
    ParametricCurve fip_mode = GetFunction(*(gloop.getTreeVector())[0], LGMIP);
    ParametricCurve fgo_mode = GetFunction(*(gloop.getTreeVector())[0], SPFGO);

    gloop.allocationAndGrowth(fip_mode,fgo_mode);
    // Command line  -writeOutput exists
    //gloop.output();
    // Prune dead parts from the trees 
    gloop.prune();
    //Set radiation use efficiency in new segments as a function of shadiness
    //experienced by mother segment
    gloop.radiationUseEfficiency();
    // Evaluate stand variables also after growth
    gloop.evaluateStandVariables();
    // collectDataAfterGrowth collects data for HDF5 file. The 0th Year dimension
    // contains initial data
    gloop.collectDataAfterGrowth(year+1);
    //Save trees as xml
    CreateTreeXMLDataSet(gloop,hdf5_trees,TXMLGROUP,gloop.getWriteInterval());
  } // End of  for(year = 0; ...)
  // [GLoop]
  /// \endinternal
  /// \par After growth collect and write results
  /// \snippet{lineno} lignum-forest.cc AGrowth
  /// \internal
  // [AGrowth]
  gloop.cleanUp();
  if (gloop.getNumberOfTrees() == 0){
    cout << "NO TREES AFTER GROWTH LOOP" <<endl;
  }
  cout << "GROWTH DONE " << "NUMBER OF TREES " << gloop.getNumberOfTrees() << endl;
  ///\par Collect HDF5 data
  ///
  ///+ Year by year, tree by tree data
  TMatrix3D<double>& hdf5_data = gloop.getHDF5TreeData();
  hdf5_file.createDataSet(TREE_DATA_DATASET_NAME,hdf5_data.rows(),hdf5_data.cols(),hdf5_data.zdim(),hdf5_data);
  hdf5_file.createColumnNames(TREE_DATA_DATASET_NAME,TREE_DATA_COLUMN_ATTRIBUTE_NAME,TREE_DATA_COLUMN_NAMES);
  ///+ Aggregate stand data
  TMatrix2D<double>& hdf5_stand_data = gloop.getHDF5StandData();
  hdf5_file.createDataSet(STAND_DATA_DATASET_NAME,hdf5_stand_data.rows(),hdf5_stand_data.cols(),hdf5_stand_data);
  hdf5_file.createColumnNames(STAND_DATA_DATASET_NAME,STAND_DATA_COLUMN_ATTRIBUTE_NAME,STAND_DATA_COLUMN_NAMES);
  ///+ Aggregate center stand data
  TMatrix2D<double>& hdf5_center_stand_data = gloop.getHDF5CenterStandData();
  hdf5_file.createDataSet(CENTER_STAND_DATA_DATASET_NAME,hdf5_center_stand_data.rows(),hdf5_center_stand_data.cols(),
			  hdf5_center_stand_data);
  hdf5_file.createColumnNames(CENTER_STAND_DATA_DATASET_NAME,STAND_DATA_COLUMN_ATTRIBUTE_NAME,STAND_DATA_COLUMN_NAMES);
  ///+ Parameters used
  TMatrix2D<double> hdf5_tree_param_data = gloop.getHDF5TreeParameterData();
  hdf5_file.createDataSet(PGROUP+TREE_PARAMETER_DATASET_NAME,hdf5_tree_param_data.rows(),hdf5_tree_param_data.cols(),
			  hdf5_tree_param_data);
  hdf5_file.createColumnNames(PGROUP+TREE_PARAMETER_DATASET_NAME,TREE_PARAMETER_ATTRIBUTE_NAME,TREE_PARAMETER_NAMES);
  ///+ Functions known in a tree
  for (unsigned int i=0; i < FN_V.size();i++){ 
    TMatrix2D<double> hdf5_tree_fn_data = gloop.getHDF5TreeFunctionData(FN_V[i]);
    hdf5_file.createDataSet(TFGROUP+FNA_STR[i],hdf5_tree_fn_data.rows(),hdf5_tree_fn_data.cols(),hdf5_tree_fn_data);
    hdf5_file.createColumnNames(TFGROUP+FNA_STR[i],TREE_FN_ATTRIBUTE_NAME,TREE_FN_COLUMN_NAMES);
  }
  ///+ All functions in the run directory
  hdf5_file.createFnDataSetsFromDir("*.fun",AFGROUP,TREE_FN_ATTRIBUTE_NAME,TREE_FN_COLUMN_NAMES);
  ///+ Command line
  vector<string> c_vec;
  std::copy( argv, argv+argc,back_inserter(c_vec));
  ostringstream cline;
  copy(c_vec.begin(),c_vec.end(),ostream_iterator<string>(cline, " "));
  hdf5_file.createDataSet(COMMAND_LINE_DATASET_NAME,cline.str());
  hdf5_file.close();
  cout << "DATA SAVED AND SIMULATION DONE" <<endl;
  // [AGrowth]
  /// \endinternal

  return 0;
}
