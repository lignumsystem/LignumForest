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
  gloop.createTrees();  
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
  /// + Create HDF5 for XML trees
  /// + Create HDF5 groups for prameters and tree functions
  /// + Create HDF5 group for XML trees
  /// \internal
  /// \snippet{lineno} lignum-forest.cc HDF5INIT
  // [HDF5INIT]
  string hdf5fname;
  ParseCommandLine(argc,argv,"-hdf5", hdf5fname);
  LGMHDF5File hdf5_file(hdf5fname);
  LGMHDF5File hdf5_trees(TREEXML_PREFIX+hdf5fname);
  hdf5_file.createGroup(PGROUP);
  hdf5_file.createGroup(TFGROUP);
  hdf5_file.createGroup(AFGROUP);
  hdf5_trees.createGroup(TXMLGROUP);
  // [HDF5INIT]
  ///@}
  /// \endinternal
  ///
  for(int year = 0; year < gloop.getIterations(); year++) {
    cout << "GROWTH LOOP YEAR " << year <<endl;
    if(gloop.getNumberOfTrees() < 5) {
      cout << "Number of trees left " << gloop.getNumberOfTrees() << " Stop." << endl;
      //Do not stop abruptly, continue to the end of loop and then write data
      //HDF5 output assumes there is at least one tree left 
      continue;
    }
    ///\defgroup C_datatransfer Save data from previous year
    ///\addtogroup lignumforest
    ///@{
    ///\par Save data
    ///+ Pass growth year to L system and to *gloop*.
    ///+ Save previous year tree height for each tree.
    ///\internal
    ///\snippet{lineno} lignum-forest.cc IGLOOP
    // [IGLOOP]
    Pine::L_age = (double)year;     //This is for L-system 
    gloop.setHPrev();
    gloop.setYear(year);
    // [IGLOOP]
    ///@}
    ///\endinternal
    ///\defgroup D_radandgrowth Radiation regime and new segments
    ///\addtogroup lignumforest
    ///@{
    ///\par The steps from the current state to a new state in the forest plot 
    ///+ Update foliage in voxel space and recalculate borderforest
    ///+ Calculate radiation climate for trees
    ///+ Calculate photosynthesis and respiration
    ///+ Create new segments
    ///+ Growth allocation
    ///+ Prune dead branches from trees
    ///\internal
    ///\snippet{lineno} lignum-forest.cc NEWSEG
    // [NEWSEG]
    gloop.setVoxelSpaceAndBorderForest();
    gloop.calculateRadiation();
    gloop.increaseXi(year);
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
  ///+ Command line
  ///+ Save HDF5 files
  ///\internal
  ///\snippet{lineno} lignum-forest.cc HDF5FILECREATE
  // [HDF5FILECREATE]
  //Tree by tree data
  TMatrix3D<double>& hdf5_data = gloop.getHDF5TreeData();
  hdf5_file.createDataSet(TREE_DATA_DATASET_NAME,hdf5_data.rows(),hdf5_data.cols(),hdf5_data.zdim(),hdf5_data);
  hdf5_file.createColumnNames(TREE_DATA_DATASET_NAME,TREE_DATA_COLUMN_ATTRIBUTE_NAME,TREE_DATA_COLUMN_NAMES);
  //Stand data
  TMatrix2D<double>& hdf5_stand_data = gloop.getHDF5StandData();
  hdf5_file.createDataSet(STAND_DATA_DATASET_NAME,hdf5_stand_data.rows(),hdf5_stand_data.cols(),hdf5_stand_data);
  hdf5_file.createColumnNames(STAND_DATA_DATASET_NAME,STAND_DATA_COLUMN_ATTRIBUTE_NAME,STAND_DATA_COLUMN_NAMES);
  //Center stand data
  TMatrix2D<double>& hdf5_center_stand_data = gloop.getHDF5CenterStandData();
  hdf5_file.createDataSet(CENTER_STAND_DATA_DATASET_NAME,hdf5_center_stand_data.rows(),hdf5_center_stand_data.cols(),
			  hdf5_center_stand_data);
  hdf5_file.createColumnNames(CENTER_STAND_DATA_DATASET_NAME,STAND_DATA_COLUMN_ATTRIBUTE_NAME,STAND_DATA_COLUMN_NAMES);
  //Parameter data
  TMatrix2D<double> hdf5_tree_param_data = gloop.getHDF5TreeParameterData();
  hdf5_file.createDataSet(PGROUP+TREE_PARAMETER_DATASET_NAME,hdf5_tree_param_data.rows(),hdf5_tree_param_data.cols(),
			  hdf5_tree_param_data);
  hdf5_file.createColumnNames(PGROUP+TREE_PARAMETER_DATASET_NAME,TREE_PARAMETER_ATTRIBUTE_NAME,TREE_PARAMETER_NAMES);
  //Tree functions
  for (unsigned int i=0; i < FN_V.size();i++){ 
    TMatrix2D<double> hdf5_tree_fn_data = gloop.getHDF5TreeFunctionData(FN_V[i]);
    hdf5_file.createDataSet(TFGROUP+FNA_STR[i],hdf5_tree_fn_data.rows(),hdf5_tree_fn_data.cols(),hdf5_tree_fn_data);
    hdf5_file.createColumnNames(TFGROUP+FNA_STR[i],TREE_FN_ATTRIBUTE_NAME,TREE_FN_COLUMN_NAMES);
  }
  //All functions
  hdf5_file.createFnDataSetsFromDir("*.fun",AFGROUP,TREE_FN_ATTRIBUTE_NAME,TREE_FN_COLUMN_NAMES);
  vector<string> c_vec;
  std::copy( argv, argv+argc,back_inserter(c_vec));
  ostringstream cline;
  copy(c_vec.begin(),c_vec.end(),ostream_iterator<string>(cline, " "));
  hdf5_file.createDataSet(COMMAND_LINE_DATASET_NAME,cline.str());
  hdf5_file.close();
  // [HDF5FILECREATE]
  ///@}
  ///\endinternal
  cout << "DATA SAVED AND SIMULATION DONE" <<endl;
  return 0;
}
