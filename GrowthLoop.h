#ifndef GROWTH_LOOP_H
/// \file GrowthLoop.h
/// \brief Implements the growth loop for LignumForest
#define GROWTH_LOOP_H
#include <cmath>
#include <cstdio>
#include <fstream>
#include <utility>
#include <mathsym.h>
#include <Point.h>
#include <PositionVector.h>
#include <TMatrix3D.h>
#include <Bisection.h>
#include <LGMHDF5File.h>
#include <Sensitivity.h> 
#include <VoxelSpace.h>
#include <TreeLocations.h>
#include <SelfThinning.h>
#include <HarvestStand.h>
#include <SomeFunctors.h>
#include <CalculateLight.h>
#include <VoxelSpace.h>
#include <DiameterGrowth.h>
#include <FindForestBoundingBox.h>
#include <XMLTree.h>
#include <ForestTrees.h>
#include <StandDescriptor.h>
#include <VoxelSpace.h>
#include <BorderForest.h>
#include <Space.h>
#include <Palubicki_functors.h>
#include <TreeDataAfterGrowth.h>
///\class CollectP
///\brief Collect photosynthates from the tree
///Use with Accumulate algorithm in stl-lignum
/// \tparam TS TreeSegment
/// \tparam BUD Bud
template <class TS, class BUD>
    class CollectP
    { 
    public:
      LGMdouble& operator()(LGMdouble &sum, TreeCompartment<TS,BUD>* tc)const
  {
    if(TS *ts = dynamic_cast<TS *>(tc))
      {
	if(GetValue(*ts,LGAWf) > 0.0) {
        sum += GetValue(*ts, LGAP);
      }
      }
    return sum;
  }
    };

/// \brief Implement the LignumForest growth loop
/// \tparam TREE Tree
/// \tparam TS Tree segment
/// \tparam BUD Bud
/// \tparam LSYSTEM Lindenmayer system
template <class TREE, class TS, class BUD,class LSYSTEM>
class GrowthLoop{
public:
  GrowthLoop()
    :vs(NULL),verbose(false),iterations(0),start_voxel_calculation(0),
     num_parts(1.0), interval(0),init_density(0),nsegment(0.0),
     tree_distance(0.0),hw_start(15),light_method(0.0),
     ws2(0.0),ws3(0.0),ws4(0.0),wr1(0.0),wfnew(0.0),
     wfaftergrowth(0.0),wsaftergrowth(0.0),qintop1(0.0),cvol(0.0),
     axisvol(0.0),qabs(0.0),treeAf(0.0),ASeg0(0.0),wood(0.0),
     wstem(0.0),wbranches(0.0),stem_sw(0.0),treeP(0.0),treeM(0.0),
     lambda(0.0),cv(0.30),to_file(false),sensitivity_analysis(false),
     crown_limit_data(false),
     writevoxels(false),increase_xi(false),
     self_thinning(false), generate_locations(false), location_file("Treelocations.txt"),
     no_trees(0), wood_voxel(true), evaluate_border_forest(true),seg_len_var(0.0),
     pairwise_self(false), eero(false),  g_fun_varies(false), g_fun_var(0.0),
     random_branch_angle(false), ba_variation(0.0) {}
  ~GrowthLoop();
  ///Initialize:  parse command line, initialize  trees, voxel space
  ///and growth loop.
  ///\param argc Numebr of command line arguments 
  ///\param argv Command line arguments
  ///\sa timeStep afterGrowth
  void initialize(int argc, char** argv);
  ///One time step
  ///\param year Number of growth steps
  void timeStep(int year);
  ///After Growth: clean up, write output
  void afterGrowth();
  ///Print usage information
  void usage()const;
  ///Command line check
  void checkCommandLine(int argc, char** argv)const;
  ///Command line parser.
  void parseCommandLine(int argc, char** argv);
  void resolveCommandLineAttributes();
  void setVerbose(bool wordy){verbose=wordy;}
  bool isVerbose()const{return verbose;}
  void initRan3(int init)const{
    int neg_init=-abs(init);
    ran3(&neg_init);
  }
  /// \brief Create trees of the forest stand.
  ///
  /// Create trees in predefined locations. Create also L-systems, data vectors and vector of file streams.
  /// Position of a tree in a tree vector defines its position in all othe vectors
  /// \pre Tree locations have been generated
  /// \sa CreateLocations
  /// \sa locations Tree locations
  /// \sa vtree vlsystem Vector of trees and their respective L-systems
  /// \sa no_h h_prev wsapwood wfoliage wroot ws_after_senescence Data vectors
  /// \sa vdatafile Tree output file streams
  void createTrees();
  /// Resize 3D array for tree data to be collected during the simulation
  /// Resize 2D array for stand data to be collected during the simulation
  /// Data will be stored in an HDF5 file after simulation.
  /// \pre Number of simulation years and trees as well as number of data columns are known
  /// \post 3D array and 2D array filled with 0:s and can be used to enter tree data
  /// \sa createTrees  parseCommandLine
  /// \sa TREE_DATA_COLUMN_NAMES STAND_DATA_COLUMN_NAMES
  /// \sa hdf5_tree_data hdf5_stand_data hdf5_center_stand_data
  void resizeTreeDataMatrix();
  void initializeTrees();
  void initializeVoxelSpace();
  void initializeFunctions();
  void initializeGrowthLoop();
  void increaseXi(int& year);
  void setTreeLocations();
  void photosynthesis(TREE& t);
  void respiration(TREE& t);
  /// \brief Collect data before growth
  ///
  /// Collect tree data for sapwood mass, foliage mass and root mass
  /// \param t tree
  /// \param i position of the tree in the tree vector
  /// \tparam TREE Lignum tree
  /// \sa wsapwood wfoliage and wroot vectors
  /// \sa vtree Tree vector
  void collectDataBeforeGrowth(TREE& t,unsigned int i);
  /// Collect tree data after growth and update data for HDF5 file
  /// \note Intial data is collected before growth loop. Thus the method should be called year=iter+1.
  /// \param year Simulation year (i.e. iteration)
  /// \pre collectDataBeforeGrowth StandDescriptor::evaluateStandVariables
  /// \sa  vtree hdf5_tree_data  hdf5_stand_data hdf5_center_stand_data wsapwood wfoliage wroot ws_after_senescence
  void collectDataAfterGrowth(const int year);
  void treeAging(TREE& t);
  double collectSapwoodMass(TREE& t);
  void setSapwoodDemandAtJunction(TREE& t);
  /// \brief Allocation of photosynthates to growth.
  /// \param t Lignum tree
  /// \param verbose Verbose output
  bool allocation(TREE& t,bool verbose);
  /// \brief Output of simulation to files.
  ///
  /// Write stand level data, target tree data, crown limit data, Fip data and the target tree xml file.
  /// \sa writeOutput writeCrownLimitData writeTreeToXMLFile writeFip
  /// \todo HDF5 implementation is advancing. Remove this method when enough data collected. Consult and
  /// agree with Risto.
  /// \sa collectDataAfterGrowth
  void output();
  /// \brief Tree level output.
  ///
  /// Write tree level output to its file. 
  /// \param t The tree
  /// \param tree_n Tree position in the tree vector
  /// \param iter Current iteration year in the similation
  /// \sa vdatafile Vector for output files for each tree
  void writeOutput(TREE& t,unsigned int tree_n,int iter);
  void writeSensitivityAnalysisData(TREE& t);
  /// \brief Write crown limit data to study crown rise.
  /// \param t The tree
  /// \param iteration Current iteration year in the similation
  void writeCrownLimitData(TREE& t,int iteration);
  /// \brief Write voxel space content (voxels)
  /// \param t The tree 
  void writeVoxels(TREE& t);
  /// \brief Write tree to file.
  /// \param t The tree
  /// \param age Tree age
  /// \param interval The write interval
  /// \pre age mod interval = 0
  void writeTreeToXMLFile(TREE& t,int age,int interval)const;
  void writeBranchInformation(TREE& t,const string& file)const;
  void writeProductionBalance(TREE& t,const string& file)const;
  /// \brief Vertical distribution of fip
  /// \param t The tree
  /// \param interval The write interval
  /// \pre Tree age mod interval = 0 and !fipfile.empty()
  /// \sa fipfile
  void writeFip(TREE& t,int interval)const;
  void prune();
  void printVariables()const;
  void cleanUp();
  void printSegmentQin();
  void printBranchMeans()const;
  void printTreeLocations(int iter)const;
  void printVoxelObjectLocations(const string& file)const;
  TREE& getTargetTree()const{
    if(target_tree > vtree.size()-1) {
      cout << "Target tree specification " << target_tree << " is wrong.";
      exit(0);}
    else
      return *(vtree[target_tree]);}
  int getTargetTreePosition()const{return (int)target_tree;}
  int getIterations()const{return iterations;}
  int getIterations() {return iterations;}
  TMatrix3D<double>& getHDF5TreeData(){return hdf5_tree_data;}
  TMatrix2D<double>& getHDF5StandData(){return hdf5_stand_data;}
  TMatrix2D<double>& getHDF5CenterStandData(){return hdf5_center_stand_data;}
  void setVoxelSpaceAndBorderForest();
  void calculateRadiation();
  StandDescriptor<TREE>& getStand() {return stand;}
  vector<TREE*>& getTrees() {return  vtree;}
  void increaseXi();
  ///\brief Photosynthesis, respiration, tree aging and data collection.
  void photosynthesisAndRespiration();
  void evaluateStandVariables();
  void createNewSegments();
  /// \brief Allocation of net photosynthates to growth.
  ///
  /// Also responsible of keeping record of tree death.
  /// Dead trees are removed from the list of trees, L-systems
  /// and from the vectors of collected data.
  /// \sa vtree vlsystem locations
  /// \sa wsapwood wfoliage wroot ws_after_senescence vdatafile
  void allocationAndGrowth();
  int getNumberOfTrees() {return no_trees;}
  void setYear(const int& y) {year = y;}
  void setHPrev(){
    for (unsigned int i = 0; i < (unsigned int)no_trees; i++)
      h_prev[(int)i] = GetValue(*vtree[i],LGAH);
  }
private:
  vector<TREE*> vtree; ///< Vector of trees
  vector<LSYSTEM*> vlsystem; ///< Vector of L-systems, one for each tree
  vector<pair<double,double> > locations; ///< Positions of trees
  vector<ofstream*> vdatafile;///< Vector of output files (as file streams) for each tree. 
  VoxelSpace *vs; ///< The voxel space spanning the forest
  StandDescriptor<TREE> stand; ///< Class to handle and print stand level quantities
  BorderForest border_forest;  ///< Homogeneous forest surrounding forest of individual trees.
  TMatrix3D<double> hdf5_tree_data; ///< 3D array[years][ntrees][ndata_cols] for trees in the stand,
                            ///< dimensions will be known after trees are generated.
  TMatrix2D<double> hdf5_stand_data; ///< 2D array[years][ndata_cols] for stand level data \sa stand
  TMatrix2D<double> hdf5_center_stand_data; ///< 2D array[years][ndata_cols] for center stand level data \sa center_stand
  bool verbose;
  bool bracket_verbose;
  int iterations; ///< Number of years simulation goes on.
  int start_voxel_calculation; ///< NOT IN USE CURRENTLY!  
  ///
  /// \brief Division of shoot foliage in dumping it to VoxelSpace
  ///
  /// Shoot is split in num_parts and foliage of each part
  /// is dumped to the corresponding voxel box. \sa setVoxelSpaceAndBorderForest
  int num_parts;
  int interval; ///< Write interval for output
  int init_density;///< NOT IN USE CURRENTLY! Init density for forced self-thinning.
  int nsegment;///< Number of TreeSegments for output
  double tree_distance ;///< Minimum allowed tree distance \sa setTreeLocations
  int  hw_start;///< Starting year of heartwood build up
  double light_method;///< NOT IN USE CURRENTLY! Mandatory method for radiation assessment
  vector<double> wsapwood; ///< Sapwood mass at start of time step 
  vector<double> wfoliage; ///< Foliage mass at start of time step
  vector<double> wroot; ///< Root mass at start of time step
  vector<double> ws_after_senescence; ///< Sapwood mass AFTER sapwood senescence
  double ws2; ///< Old segment sapwood mass
  double ws3; ///< Sapwood mass allocated to diameter growth
  double ws4; ///< Sapwood mass in the new segments
  double wr1; ///< Root mass required by new foliage
  double wfnew; ///<  New foliage mass
  double wfaftergrowth;///< Foliage after new growth
  double wsaftergrowth;//< Sapwood after new growth
  double qintop1; ///< Incoming radioation at tree top
  double cvol; ///< Crown volume of one tree \sa cv
  double axisvol;///< Volume of main stem
  double qabs; ///< Intercepted radiation of a tree
  double treeAf; ///< Foliage area of the tree
  double ASeg0; ///< Surface area of new segments
  double wood; ///< Wood mass of the tree
  double wstem; ///< Wood mass in the main axis (stem)
  double wbranches; ///< Wood mass in branches
  double stem_sw; ///< Sapwood mass in the main stem
  ///
  /// \brief Photosynthetic production of the tree.
  ///
  ///Careful when  using this  variable.  After the vector  of trees  has been
  ///passed over, treeP holds the value
  ///of the last tree in the vector \sa allocationAndGrowth
  double treeP;
  double treeM; ///< Respiration of tree. \warning The same warning as for treeP \sa TreeP

  ///
  /// \brief lambda = Iteration parameter of new growth
  ///
  /// After allocation Photosynthesis - Respiration - Growth(lambda) = 0
  /// \sa allocationAndGrowth
  double lambda;
  summing bs;///< Mean branch length
  DCLData dcl;///< Diameter and heigth at the crown base.

  /// \brief This functor class returns crown volume of a tree
  ///
  /// Class CrownVolume declared and defined in stl-lignum/TreeFunctor.h
  /// and stl-lignum/TreeFunctorI.h.
  CrownVolume<TS,BUD> cv;

  /// \brief This functor class returns stem volume of a tree
  ///
  /// Class MainAxisVolume Defined in stl-lignum/TreeFunctor.h
  MainAxisVolume<TS,BUD> mav;
  string metafile; ///< File of names of actual parameter files etc. \sa checkCommandLine
  string voxelfile;///< File of parameters for the VoxelSpace \sa checkCommandLine
  string xmlfile; ///< XML file where the tree can be saved and restored from \sa parseCommandline
  string datafile;///< Text file to save tree data each time step \sa parseCommandline

  /// \brief File name for storing radiation condition information of shoots
  ///
  /// The radiation conditions of a shoot are determined with the aid of relative
  /// incoming radiation = incoming radiation to shoot / max.  incoming radiation, and
  /// further, in terms of fip = f(relative incoming radiation). The fip function affects
  /// length of new shoot and has been defined in Perttunen et al. 1996. \sa writeFip
  string fipfile;
  bool to_file;///< If output to a file \sa createTrees   sa\ writeOutput
  bool sensitivity_analysis; ///< If sensitivity analysis is done. \sa parseCommandLine
  bool crown_limit_data; ///< If output about crown base. \sa parseCommandLine \sa writeCrownLimitData
  bool writevoxels; ///< If output about the voxelspace \sa writeVoxels

  /// \brief If the primary wood proportion (denoted xi) in new shoots increases
  ///
  /// The effect of primary wood for sapwood proportion in new segments is described in
  /// Perttunen et al. 1996. Primary wood proportion = LGPXi LIGNUM parameter, see
  /// stl-lignum/include/LGMSymbols.h \sa increaseXi
  bool increase_xi;
  int xi_start; ///< Starting year of Xi increase \sa increase_xi
  bool self_thinning;///< NOT IN USE CURRENTLY! If forced self thinning.

  /// \brief The function K appears in radiation interception calculation.
  ///
  /// The K function is defined in Perttunen et al. 1998, Eq. 8
  /// It is in action in stl-lignum/include/ShadingI.h
  /// \sa EvaluateRadiationForCfTreeSegmentInVoxelSpace
  ParametricCurve K;
  ParametricCurve stems_ha;///< Density of the forest as a function of age if it is forced
  ParametricCurve fdensity_ha;///< Density of border forest as a f. of age if forced
  Sensitivity<TS,BUD> sensitivity; ///< For printing out sensitivity analysis results
  bool generate_locations; ///< If tree positions are set by the program. \sa setTreeLocations
  string location_file; ///< Tree positions are read from here. \sa setTreeLocation
  ifstream location_stream; ///< SEEMS SUPERFLUOUS!
  int no_trees; ///< Number of trees in the forest
  bool noWoodVoxel; ///< SUPERFLUOUS!
  bool wood_voxel; ///< If woody parts are dumped into voxels \sa setVoxelSpaceAndBorderForest
  int year; ///< Year in simulation (start = 0)
  unsigned int target_tree; ///< Information of this tree is printed out. \sa output
  bool write_output; ///< If information of target tree is printed out. \sa output
  ofstream* stand_output; ///< Strem for target tree output \sa output
  StandDescriptor<TREE> center_stand; ///< To deal with center part of stand \sa setTreeLocations
  ofstream* cstand_output; ///< Stream for center stand output
  bool evaluate_border_forest; ///< If border forest in radiation calculations? \sa calculateRadiation
  LGMdouble k_border_conifer;///< Extinction coeffient for border forest conifers \sa calculateRadiation
  /// \brief To keep track of trees that do not grow in height
  ///
  /// If a tree does not grow for a while it is considered dead.
  /// This variable keeps track of non-growing trees. \sa h_prev \sa allocationAndGrowth
  vector<int> no_h;
  vector<double> h_prev; ///< Tree height of previous time step. \sa no_h \sa allocationAndGrowth
  LGMdouble p0_var; ///< Random variation in photosynthetic efficiency between trees. \sa initializeTree
  string stand_file; ///< File name for output of stand values \sa output
  string cstand_file; ///< File name for output of center stand values \sa output
  /// \brief For generation of random variation in lengths of new segments.
  ///
  /// The variability is controlled with LIGNUM parameter LGPlen_random (see
  /// stl-lignum/include/LGMSymbols.h) \sa SetScotsPineSegmentLength
  LGMdouble seg_len_var; 
  //===============  24.10.2019
  /// \brief Shading in tree's own crown is analyzed with backward ray casting.
  ///
  /// Backward ray casting (leads to pairwise comparison of shoots) is explained in
  /// Perttunen et al. 1998 and Sievanen et al. 2008. \sa calculateRadiation
  bool pairwise_self;
  /// \brief For studying the effects of variation in tree's physiological properties.
  ///
  /// Variation in parameter values makes it possible to study the growth effects of
  /// a) leaf light climate and b) syncronized variation in leaf nitrogen concentration,
  /// leaf mass per area and leaf longevity, shoot length, and shoot leaf area and,
  /// axis thickness scaling from the base of the stem to the axis tips. (-eero
  /// stems from Eero Nikinmaa who proposed this) \sa initializeTrees
  bool eero;

  /// \brief Variability in the function: shoot growth = f(branching order)
  ///
  /// Shoot growth normally depends on the branching order of the growing shoot.
  /// This function is explained in Sievanen et al. 2018 Eq. 4. This option
  /// varies the function randomly between trees. \sa initializeTrees
  bool g_fun_varies;
  double g_fun_var;  ///< Amount of random variation in f(branching order) \sa g_fun_varies
  bool random_branch_angle; ///< Branching angle varies in trees. \sa initializeTrees
  double ba_variation;  ///< Amount of variation in braching angle. \sa random_branch_angle
  //=================== 30.9.2021
  /// \brief Shoot growth is controlled by Extended Borchert-Honda mechanism.
  ///
  ///The mechanism is explained in Sievanen et al. 2018 Eqs 6 and 7.
  /// \sa createNewSegments \sa SetScotsPineSegmentLength
  bool growthloop_is_EBH;
  bool growthloop_is_EBH1; ///< As growthloop_is_EBH but with one parameter. \sa  growthloop_is_EBH
  double growthloop_EBH1_value; ///< Parameter value for growthloop_is_EBH1 \sa  growthloop_is_EBH1
};
#endif
#include <GrowthLoopI.h>
