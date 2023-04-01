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
using namespace LignumForest;
using namespace PineTree;
///\defgroup mainglobals Global variables for growth allocation
///Declaration of a number of global variables for (*ad hoc*) growth allocation.
///They are easy to add to the program but maybe should be implemented as function arguments or something else,
///see e.g. PineBudData in Pine.h \sa PineBudData
/// @{
///If function adhoc is in use: it boosts shoot growth in lower parts of the crown, see Eq. 8 in Sievanen et al. 2018 \sa adhoc
bool is_adhoc = false;
///This function boosts shoot growth in lower parts of the crown, see Eq. 8 in Sievanen et al. 2018
///\sa SetScotsPineSegmentLength
ParametricCurve adhoc("adhoc.fun");
///Height of grown base for function adhoc \sa adhoc
double global_hcb=0.0;
///@}
///\defgroup spaceo Global variables to asses space occupation in shoot growth
/// @{
///Firmament for calculation if foliage in growth direction (Eq 10 in Sievanen et al. 2018).
///\sa SetScotsPineSegmentLength
Firmament dummy_firm;
///VoxelSpace for calculation if foliage in growth direction
///Eq 10 in Sievanen et al. 2018 \sa SetScotsPineSegmentLength
VoxelSpace space_occupancy(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
			   0.1,0.1,0.1,5,5,5,dummy_firm);
///Only one voxelbox at the end of new Segment is checked
bool space0 = false;
///Also neighboring boxes in direction of new Segment are checked
bool space1 = false;
///Also all neighboring boxes are checked
bool space2 = false;
///Search distance of neighboring boxes for space2
double space2_distance = 0.3;
///@}
///\defgroup spaceb Global variables to asses amount of foliage in front of a bud
/// @{ 
///`bud_view_f` global function is used to asses foliage area density in the view cone of a bud. It affects
///the number of lateral buds created in L-system. Eq 11 in Sievanen et al. 2018.
///\sa pine-em98.L  \sa SetBudViewFunctor
ParametricCurve bud_view_f;
///Sets bud view function in use. \sa bud_view_f 
bool is_bud_view_function = false;   
///@}
///\defgroup Lextglob Global L-system variables
///Global variables to pass information to and from L-system
///@{
///Conveys tree age to the L-system \sa pine-em98.L \sa L_age
double L_age = 0.0;
///Conveys tree height to the L-system \sa pine-em98.L \sa L_H
double L_H =0.0;
///@}
///\defgroup random Random seed
///@{
///Random seed for `ran3` function. \sa GrowthLoop
int ran3_seed=-1;
///@}
///\defgroup  ranvar Global variables for tree variation
///Global variables to generate random variation between tree individuals.
///They work in the Lindenmayer system \sa pine-em98.L
/// @{
///Initial height ot trees
double H_0_ini=0.0;
///Variation of initial heights
double  H_var_ini=0.0;
///Variation in the initial number of buds.
int n_buds_ini_min=0.0;
int n_buds_ini_max=0.0;
///Variation in the number of buds
double rel_bud=0.0;
///If bud variation is on \sa pine-em98.L
bool bud_variation=false;
///Variation in branching_angle
double branch_angle=0.0; 
///@}

///\typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud, Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest
///ScotsPineForest (i.e. GrowthLoop template instance) captures the growth loop of a forest stand.
typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,
		   Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > ScotsPineForest;

/// \defgroup groupmain Main program for growth loop
/// @{
/// Content of the main program to run Growth loop for LignumForest.
/// \fn int main(int argc, char** argv)
/// \brief Main function for growth.
/// Check and parse command line. The  command line includes switches to control
/// the simulation.
/// \param argc An integer counter number of command line arguments.
/// \param argv Vector of command line argument strings.
int main(int argc, char** argv)
{
  /// **Growth loop variables**
  /// \snippet{lineno} lignum-forest.cc Vars
  /// \internal
  // [Vars]
  Sensitivity<ScotsPineSegment,ScotsPineBud> sensitivity;                                                           
  GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> > gloop;
  // [Vars]
  /// \endinternal
  
  ran3(&ran3_seed);

  // \snippet{lineno} lignum-forest.cc Mf
  // \internal
  // MixedForest  will  implement a  class  that  can  manage a  forest
  // consting of several tree species.
  // This is a quick way to test that the class compiles.
  // [Mf] 
  // MixedForest<ScotsPineForest,ScotsPineForest> mf;
  // mf.initialize(argc,argv);
  // mf.growthLoop();
  // mf.afterGrowth();
  // [Mf]
  /// \endinternal
  ///**Initialize forest**
  ///\snippet{lineno} lignum-forest.cc InitForest
  ///\internal
  // [InitForest]
  gloop.parseCommandLine(argc,argv);
  gloop.resolveCommandLineAttributes();
  gloop.printVariables();
  //Reads in some functions from files
  gloop.initializeFunctions();
  //Configure tree locations
  gloop.setTreeLocations();
  //Depending on configuration random or pretermined locations
  gloop.createTrees();  
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
  //Initialize  HDF5 content
  string hdf5fname;
  ParseCommandLine(argc,argv,"-hdf5", hdf5fname);
  LGMHDF5File hdf5_file(hdf5fname);
  LGMHDF5File hdf5_trees(TREEXML_PREFIX+hdf5fname);
  hdf5_file.createGroup(PGROUP);
  hdf5_file.createGroup(TFGROUP);
  hdf5_file.createGroup(AFGROUP);
  hdf5_trees.createGroup(TXMLGROUP);
  /// **Growth loop**
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
    gloop.allocationAndGrowth();
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
    //Save as xml
    CreateTreeXMLDataSet(gloop,hdf5_trees,TXMLGROUP,gloop.getWriteInterval());
  } // End of  for(year = 0; ...)
  // [GLoop]
  /// \endinternal
  /// **After growth collect and write results**
  /// \snippet{lineno} lignum-forest.cc AGrowth
  /// \internal
  // [AGrowth]
  gloop.cleanUp();
  if (gloop.getNumberOfTrees() == 0){
    cout << "NO TREES AFTER GROWTH LOOP" <<endl;
  }
  cout << "GROWTH DONE " << "NUMBER OF TREES " << gloop.getNumberOfTrees() << endl;
  //Collect HDF5 data
  //Year by year, tree by tree data
  TMatrix3D<double>& hdf5_data = gloop.getHDF5TreeData();
  hdf5_file.createDataSet(TREE_DATA_DATASET_NAME,hdf5_data.rows(),hdf5_data.cols(),hdf5_data.zdim(),hdf5_data);
  hdf5_file.createColumnNames(TREE_DATA_DATASET_NAME,TREE_DATA_COLUMN_ATTRIBUTE_NAME,TREE_DATA_COLUMN_NAMES);
  //Aggregate stand data
  TMatrix2D<double>& hdf5_stand_data = gloop.getHDF5StandData();
  hdf5_file.createDataSet(STAND_DATA_DATASET_NAME,hdf5_stand_data.rows(),hdf5_stand_data.cols(),hdf5_stand_data);
  hdf5_file.createColumnNames(STAND_DATA_DATASET_NAME,STAND_DATA_COLUMN_ATTRIBUTE_NAME,STAND_DATA_COLUMN_NAMES);
  //Aggregate center stand data
  TMatrix2D<double>& hdf5_center_stand_data = gloop.getHDF5CenterStandData();
  hdf5_file.createDataSet(CENTER_STAND_DATA_DATASET_NAME,hdf5_center_stand_data.rows(),hdf5_center_stand_data.cols(),
			  hdf5_center_stand_data);
  hdf5_file.createColumnNames(CENTER_STAND_DATA_DATASET_NAME,STAND_DATA_COLUMN_ATTRIBUTE_NAME,STAND_DATA_COLUMN_NAMES);
  //Parameters used 
  TMatrix2D<double> hdf5_tree_param_data = gloop.getHDF5TreeParameterData();
  hdf5_file.createDataSet(PGROUP+TREE_PARAMETER_DATASET_NAME,hdf5_tree_param_data.rows(),hdf5_tree_param_data.cols(),
			  hdf5_tree_param_data);
  hdf5_file.createColumnNames(PGROUP+TREE_PARAMETER_DATASET_NAME,TREE_PARAMETER_ATTRIBUTE_NAME,TREE_PARAMETER_NAMES);
  //Functions known in a tree
  for (unsigned int i=0; i < FN_V.size();i++){ 
    TMatrix2D<double> hdf5_tree_fn_data = gloop.getHDF5TreeFunctionData(FN_V[i]);
    hdf5_file.createDataSet(TFGROUP+FNA_STR[i],hdf5_tree_fn_data.rows(),hdf5_tree_fn_data.cols(),hdf5_tree_fn_data);
    hdf5_file.createColumnNames(TFGROUP+FNA_STR[i],TREE_FN_ATTRIBUTE_NAME,TREE_FN_COLUMN_NAMES);
  }
  //All functions used
  hdf5_file.createFnDataSetsFromDir("*.fun",AFGROUP,TREE_FN_ATTRIBUTE_NAME,TREE_FN_COLUMN_NAMES);
  //Command line
  vector<string> c_vec;
  std::copy( argv, argv+argc,back_inserter(c_vec));
  ostringstream cline;
  copy(c_vec.begin(),c_vec.end(),ostream_iterator<string>(cline, " "));
  hdf5_file.createDataSet(COMMAND_LINE_DATASET_NAME,cline.str());
  hdf5_file.close();
  cout << "DATA SAVED AND SIMULATION DONE" <<endl;
  //gloop.writeTreeToXMLFile(gloop.getTargetTree(),GetValue(gloop.getTargetTree(),LGAage),1);
  // [AGrowth]
  /// \endinternal

// #if defined (__APPLE__) || defined(__MACOSX__)
//   if(CheckCommandLine(argc,argv,"-viz")) {
//     LGMVisualization viz;
//     viz.InitVisualization(argc,argv);
//     // textures 512x512
//     viz.AddCfTree(gloop.getTargetTree(), "Manty.bmp", "neulaset5.tga");
//     float th = (float)GetValue(gloop.getTargetTree(),LGAH);
//     cout << th << endl;
//     //viz.ResetCameraPosition(th);
//     viz.SetMode(SOLID);
//     //viz.ResetCameraPosition(GetValue(gloop.getTargetTree(),LGAH));
//     viz.StartVisualization();
//   }
// #endif
  return 0;
}
/// @} Group main
