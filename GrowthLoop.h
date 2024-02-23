#ifndef GROWTH_LOOP_H
/// \file GrowthLoop.h
/// \brief Implements the growth loop for LignumForest
#define GROWTH_LOOP_H
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
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
#include <DiameterGrowth.h>
#include <FindForestBoundingBox.h>
#include <XMLTree.h>
#include <ForestTrees.h>
#include <StandDescriptor.h>
#include <BorderForest.h>
#include <Space.h>
#include <Palubicki_functors.h>
#include <TreeDataAfterGrowth.h>

namespace LignumForest{
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
  class GrowthLoop:public LGMHDF5{
    template <class TREE1,class TS1,class BUD1,class LSYSTEM1>
    /// \brief Create XML string reprentations for the trees in the forest stand.
    ///
    /// Each tree will be its own dataset in its age group. Trees are (must be) collected
    /// during the growth loop.
    /// \pre The HDF5 file `hdf5_file` must be open
    /// \param gl The GrowLoop
    /// \param hdf5_file The HDF5 file where the XML strings will be stored
    /// \param group_name The name of the root group (dataset) for the XML strings
    /// \param interval Write interval: trees will be written when `age mod interval = 0`.
    /// \return 0 Always returns zero. 
    /// \note The dataset naming for the trees will be */dataset_name/`age`/Tree_`tree_id`*,
    /// where the `age` is the age of the tree and `tree_id` the ID of the tree.
    /// \post  The HDF5 file `hdf5_file` remains open
    /// \attention The `hdf5_file` must be closed after the growth loop before program exit.
    /// \todo Improve the use of the return value to use return values from HDF5 functions
    friend int CreateTreeXMLDataSet(const GrowthLoop<TREE,TS,BUD,LSYSTEM>& gl, LGMHDF5File& hdf5_file,const string& group_name,
				    const int interval);
    /// \brief Update `ws_after_senescence` vector
    /// \param gl GrowthLoop
    /// \param ws Sapwood mass
    /// \param pos position in the `ws_after_senescence` vector
    /// \sa ws_after_senescence vector
    friend void UpdateSapwoodAfterSenescence(GrowthLoop<TREE,TS,BUD,LSYSTEM>& gl,double ws,unsigned int pos)
    {
      gl.ws_after_senescence[pos]=ws;
    }
  public:
    ///\brief Initialize variables.
    ///
    ///No tree. The `location_file` is set to "TreeLocations.txt" for predefined tree locations.   
    ///`hw_start` is set to 15.   
    ///`cv` is set to 0.30.      
    ///\sa hw_start cv
    GrowthLoop()
      :vs(NULL),verbose(false),iterations(0),start_voxel_calculation(0),
       num_parts(1.0), interval(0),init_density(0),nsegment(0.0),
       tree_distance(0.0),hw_start(15),light_method(0.0),
       ws2(0.0),ws3(0.0),ws4(0.0),wr1(0.0),wfnew(0.0),
       wfaftergrowth(0.0),wsaftergrowth(0.0),qintop1(0.0),cvol(0.0),
       axisvol(0.0),qabs(0.0),treeAf(0.0),ASeg0(0.0),wood(0.0),
       wstem(0.0),wbranches(0.0),stem_sw(0.0),treeP(0.0),treeM(0.0),
       cv(0.30),to_file(false),sensitivity_analysis(false),
       crown_limit_data(false),
       writevoxels(false),increase_xi(false),
       self_thinning(false), generate_locations(false), location_file("Treelocations.txt"),
       no_trees(0), stand_output(NULL),cstand_output(NULL),wood_voxel(true), evaluate_border_forest(true),seg_len_var(0.0),
       pairwise_self(false), eero(false),  g_fun_varies(false), g_fun_var(0.0),
       random_branch_angle(false), ba_variation(0.0), dDb(0.003) {}
    ///\brief Initialize with one tree
    ///
    ///One tree with its L_system. The `location_file` is set to "TreeLocations.txt" for predefined tree locations.   
    ///`hw_start` is set to 15.       
    ///`cv` is set to 0.30.    
    ///`no_tree is set 1 .  
    ///`wsapwood`, `wfoliage`, `wroot`, ws_after_senescence` and `lambdav` vectors are length 1 inital value 0.0.
    ///\sa hw_start cv
    ///\attention Meant only to collect HDF5 tree data. No growth.
    ///\remark Used in CrownDensity project to collect HDF5 data with existing methods.
    GrowthLoop(TREE* t, LSYSTEM* l)
      :vs(NULL),verbose(false),iterations(0),start_voxel_calculation(0),
       num_parts(1.0), interval(0),init_density(0),nsegment(0.0),
       tree_distance(0.0),hw_start(15),light_method(0.0),
       ws2(0.0),ws3(0.0),ws4(0.0),wr1(0.0),wfnew(0.0),
       wfaftergrowth(0.0),wsaftergrowth(0.0),qintop1(0.0),cvol(0.0),
       axisvol(0.0),qabs(0.0),treeAf(0.0),ASeg0(0.0),wood(0.0),
       wstem(0.0),wbranches(0.0),stem_sw(0.0),treeP(0.0),treeM(0.0),
       cv(0.30),to_file(false),sensitivity_analysis(false),
       crown_limit_data(false),
       writevoxels(false),increase_xi(false),
       self_thinning(false), generate_locations(false), location_file("Treelocations.txt"),
       no_trees(1),stand_output(NULL),cstand_output(NULL),wood_voxel(true), evaluate_border_forest(true),seg_len_var(0.0),
       pairwise_self(false), eero(false),  g_fun_varies(false), g_fun_var(0.0),
       random_branch_angle(false), ba_variation(0.0), dDb(0.003)
    {
      vtree.push_back(t);
      vlsystem.push_back(l);
      wsapwood.push_back(0.0); 
      wfoliage.push_back(0.0); 
      wroot.push_back(0.0); 
      ws_after_senescence.push_back(0.0);
      lambdav.push_back(0.0);
    }
    ~GrowthLoop();
    ///\brief Initialize based on command line
    ///
    /// Parse command line, initialize trees and functions, voxel space and growth loop.
    /// \attention Not implemented
    ///\param argc Numebr of command line arguments 
    ///\param argv Command line arguments
    /// \todo Collect parsing command line, initializing functions, creating and initializing trees and voxel space
    /// as well as initialising growth loop to this method.
    ///\sa The sequence of parseCommandLine resolveCommandLineAttributes initializeFunctions
    ///\sa setTreeLocations createTrees initializeTrees initalizeVoxelSpace resizeTreeDataMatrix
    ///\sa initializeGrowthLoop in main loop      
    void initialize(int argc, char** argv);
    ///One time step
    ///\param year Number of growth steps
    void timeStep(int year);
    ///After Growth: clean up, write output
    void afterGrowth();
    ///Print usage information
    ///\deprecated Options for output files 
    ///\remark Simulation results as well as XML trees will be in HDF5 files. See option `-hdf5` 
    void usage()const;
    ///\brief Command line check
    void checkCommandLine(int argc, char** argv)const;
    ///\brief Command line parser.
    void parseCommandLine(int argc, char** argv);
    ///\brief Currently checks `eero` from command line
    ///\post If `eero` = *true* then `bud_variation` = *false*
    ///\sa eero
    void resolveCommandLineAttributes();
    ///\brief verbose output
    ///\param wordy If *true* then verbose output
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
    /// \sa setTreeLocations
    /// \sa locations Tree positions
    /// \sa vtree vlsystem Vector of trees and their respective L-systems
    /// \sa no_h h_prev wsapwood wfoliage wroot ws_after_senescence Data vectors
    /// \sa vdatafile Tree output file streams
    /// \note The command line option `-generateLocations <num>` will override all other forest generation options 
    void createTrees();
    /// \brief Resize data structures for HDF5 file.
    ///
    ///  Resize data structures for HDF5 file to their right sizes before simulation.
    /// -# Resize TMatrix 3D array for tree data to be collected during the simulation.
    /// -# Resize TMatrix 2D arrays for stand and center stand data to be collected during the simulation.    
    ///
    /// Data will be stored in an HDF5 file after simulation.
    /// \pre Simulation years and number of trees at the beginning as well as number of data columns are known
    /// \post 3D array and 2D array filled with std:nan and can be used to enter tree data
    /// \sa createTrees  parseCommandLine
    /// \sa TREE_DATA_COLUMN_NAMES STAND_DATA_COLUMN_NAMES
    /// \sa hdf5_tree_data hdf5_stand_data hdf5_center_stand_data
    void resizeTreeDataMatrix();
    void initializeTrees();
    void initializeVoxelSpace();
    /// \brief Read ParametricCurve files
    /// \pre ParametricCurve files must exists
    /// \todo Divide ParametricCurve files to mandatory (throw exception) and optional (warning to standard error)
    /// \todo Replace hard-coded file names with parameterised
    /// \remark If `verbose` output the health of parametric curves is echoed with ParametricCurve::ok
    /// \remark The existence of a file is tested with `ìfstream`. In the C++ 17 STL
    /// one can clean up implementation with std::filesystem::exists
    void initializeFunctions();
    /// \brief Initialize L-systems and Tree root masses
    void initializeGrowthLoop();
    /// \brief Increase the value of Lignum::LGPxi in trees after given year
    /// \param year The start year to increase Lignum::LGPxi
    /// \attention Currently hard coded without specific parametrization
    /// \note The implementation is the same as in CrownDensity
    void increaseXi(int& year);
    ///Generate random tree locations, or read them from a file
    ///Establish also stand corners with this information
    ///They are set in StandDescriptor and BorderForest
    ///\sa GenerateLocations
    void setTreeLocations();
    void photosynthesis(TREE& t);
    void respiration(TREE& t);
    /// \brief Collect data before growth.
    ///
    /// Collect tree data for sapwood mass, foliage mass and root mass
    /// \tparam TREE Lignum tree
    /// \param t the tree
    /// \param i position of the tree in the tree vector
    /// \sa wsapwood wfoliage and wroot vectors
    /// \sa vtree Tree vector
    void collectDataBeforeGrowth(TREE& t,unsigned int i);
    /// \brief Collect tree data after growth and update data for HDF5 file.
    ///
    /// Collect data for each tree for a single year to `tdafter` dictionary (STL data type *map*).
    /// Add a row for each tree in 3D matrix `hdf5_tree_data`.  
    /// Collect forest stand and center stand level aggregate data into 2D arrays `sdafter` and `csdafter` respectively from the forest stand.  
    /// Intial data is collected before growth loop filling the first (0th) year.
    /// In this case the method should be called parameter `year`=iter+1.
    /// \param year Simulation year (i.e. iteration)
    /// \param collect_stand If *true* collect forest stand and center stand level aggregate data
    /// \post Each tree maintains its position (row) in 3D hdf5_tree_data denoted by TreeId number.
    /// \sa GrowthLoop::collectDataBeforeGrowth StandDescriptor::evaluateStandVariables GrowthLoop::resizeTreeDataMatrix
    /// \sa hdf5_tree_data  hdf5_stand_data hdf5_center_stand_data Data for HDF5 file
    /// \sa wsapwood wfoliage wroot ws_after_senescence Collect biomasses  
    /// \sa vtree Tree vector 
    void collectDataAfterGrowth(const int year, bool collect_stand=true);
    /// \brief Collect sapwood after senescence.
    ///
    /// Collect sapwood mass to `ws_after_senescence` vector.
    /// \pre treeAging must have been called 
    /// \param t The Lignum tree
    /// \param i position of the tree in `vtree` vector
    /// \sa collectSapwoodMass
    void collectSapwoodAfterSenescence(TREE& t, unsigned int i);
    /// TreeAging takes care of senescence in segments (above ground part) and root mortality.
    /// \param t Lignum tree
    void treeAging(TREE& t);
    double collectSapwoodMass(TREE& t);
    void setSapwoodDemandAtJunction(TREE& t);
    /// Allocation of photosynthates to growth.
    /// The Allocation of photosynthates P after respiration costs M,  that is,
    /// finding value of lambda parameter that makes demand, P-M, to match available
    /// resources G with the aid of iteration:
    ///    P - M = G(lambda)
    /// that can be use in growth.
    ///
    /// \return true if lambda s.t. P-M-G(lambda)=0, false if  P - M < 0 or iteration cannot find solution
    /// \remark If false the tree `t` is considered dead and will be removed from tree vector 
    /// \tparam TREE Lignum tree
    /// \param t Lignum tree 
    /// \param verbose Verbose output
    /// \param fip_mode Function fip(ip) after growth mode change
    /// \param fgo_mode Function fgo(go) after growth mode change
    /// \sa vtree
    bool allocation(TREE& t,bool verbose,const ParametricCurve& fip_mode, const ParametricCurve& fgo_mode);
    /// \brief Set radiation use efficiency (rue) in new segments
    /// Set the radiation use efficiency (rue) of new segments (age = 0) on the basis of shadiness
    /// experienced by their mother. It is measured as Qin / QinMax, QinMax = maximum Qin
    /// in the forest = diffuseBallSensor() reading of Firmament (it is the same for all trees in the forest).
    void radiationUseEfficiency();
    /// \brief Output of simulation to files.
    ///
    /// Write stand level data, target tree data, crown limit data, Fip data and the target tree xml file.
    /// \sa writeOutput writeCrownLimitData writeTreeToXMLFile writeFip
    /// \deprecated HDF5 implementation is advancing. Remove this method when enough data collected. Consult and
    /// agree with Risto.
    /// \sa collectDataAfterGrowth
    void output();
    /// \brief Tree level output.
    ///
    /// \deprecated Write now tree level output to HDF5 file. 
    /// \tparam TREE Lignum tree
    /// \param t The tree
    /// \param tree_n Tree position in the tree vector
    /// \param iter Current iteration year in the similation
    /// \sa vdatafile Vector for output files for each tree
    void writeOutput(TREE& t,unsigned int tree_n,int iter);
    void writeSensitivityAnalysisData(TREE& t);
    /// \brief Write crown limit data to study crown rise.
    /// \tparam TREE Lignum tree
    /// \param  t The tree 
    /// \param iteration Current iteration year in the similation
    void writeCrownLimitData(TREE& t,int iteration);
    /// \brief Write voxel space content (voxels)
    /// \tparam t The tree 
    void writeVoxels(TREE& t);
    /// \brief Write tree to file.
    /// \tparam TREE Lignum tree
    /// \param t The tree
    /// \param age Tree age
    /// \param interval The write interval
    /// \pre age mod interval = 0
    void writeTreeToXMLFile(TREE& t,int age,int interval)const;
    void writeBranchInformation(TREE& t,const string& file)const;
    void writeProductionBalance(TREE& t,const string& file)const;
    /// \brief Vertical distribution of fip
    /// \tparam TREE Lignum  tree
    /// \param t The tree
    /// \param interval The write interval
    /// \pre Tree age mod interval = 0 and !fipfile.empty()
    /// \sa fipfile
    void writeFip(TREE& t,int interval)const;
    void prune();
    void printVariables()const;
    /// \brief Collect data if èero` = *true* to *eero.dat*
    ///
    /// Collect data to *eero.dat*. Call L-system End function. Close stand and center stand output files.
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
    int getWriteInterval()const{return interval;}
    ///3D Data array[year][tree][data_cols] for each year for each tree 
    TMatrix3D<double>& getHDF5TreeData(){return hdf5_tree_data;}
    ///2D data array[year][data_cols] for stand level aggregate data
    TMatrix2D<double>& getHDF5StandData(){return hdf5_stand_data;}
    ///2D data array[year][data_cols] aggregate data for the center part of the stand
    TMatrix2D<double>& getHDF5CenterStandData(){return hdf5_center_stand_data;}
    ///Tree parameter values used in simulation
    TMatrix2D<double> getHDF5TreeParameterData();
    ///Find function definition known to a tree 
    ///\pre Tree vector must have trees
    ///\return TMatrix2D<double>(N,2) of (x,f(x)) values (N rows, 2 columns)
    ///\sa cxxadt::ParametricCurve LGMF vtree
    TMatrix2D<double> getHDF5TreeFunctionData(const LGMF fn_enum);
    void setVoxelSpaceAndBorderForest();
    ///Calculates radiation for the trees 
    ///-# If *pairwise_self* pairwise for a tree itself
    ///-# Through voxelspace 
    ///-# Account for border forest
    ///\sa pairwise_self
    ///\sa Lignum::EvaluateRadiationForCfTreeSegmentInVoxelSpace
    void calculateRadiation();
    StandDescriptor<TREE>& getStand() {return stand;}
    ///\brief Access the vector of trees.
    ///
    ///\return `vtree` the vector of trees
    ///\sa vtree
    const vector<TREE*>& getTrees() const  {return  vtree;}
    ///\brief Photosynthesis, respiration, tree aging and data collection.
    ///
    /// Photosynthesis, respiration, tree aging and data collection for all trees in `vtree`.
    /// \attention The following steps are taken:
    /// -# photosynthesis
    /// -# respiration
    /// -# collectDataBeforeGrowth
    /// -# treeAging
    /// -# collectSapwoodMass
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
    void allocationAndGrowth(const ParametricCurve& fip_mode, const ParametricCurve& fgo_mode);
    int getNumberOfTrees() {return no_trees;}
    void setYear(const int& y) {year = y;}
    void setHPrev(){
      for (unsigned int i = 0; i < (unsigned int)no_trees; i++)
	h_prev[(int)i] = GetValue(*vtree[i],LGAH);
    }
    /// \brief Harvest: remove percentage of trees from the forest stand 
    /// \param percentage Percentage of shortest trees (number of trees) to be removed
    /// \post The tree vector `vtree` and associated data vectors updated
    /// \note Other criteria like several harvest times and harvesting based on basal area may follow
    /// \sa removeTreesAllOver
    void harvestForest(double percentage);
    /// \brief Remove trees from the forest stand.
    /// \param vremove Positions of trees in `vtree` to be removed.
    /// \pre Tree positions in `vremove` must be in ascending order
    /// \post The tree vector `vtree` and associated data vectors updated
    /// \sa harvestForest
    /// \sa vtree vlsystem locations
    /// \sa wsapwood wfoliage wroot ws_after_senescence vdatafile
    void removeTreesAllOver(const vector<unsigned int>& vremove);
    vector<TREE*>& getTreeVector() {return vtree;}
  private:
    vector<TREE*> vtree; ///< Vector of trees. \sa getTrees
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
    ///
    /// \brief lambda = Iteration parameter of new growth
    /// Save lambda for each tree.
    /// After allocation Photosynthesis - Respiration = Growth(lambda)
    /// \sa allocationAndGrowth
    vector<double>  lambdav;
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
    int init_density;///< Init density for forced self-thinning. \attention Not in use currently! 
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
    /// stl-lignum/include/LGMSymbols.h \sa xi_start
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
    ParametricCurve fip_mode;///< Function fip after growth mode change
    ParametricCurve fgo_mode;///< Function fgo after growth mode change
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
    bool growthloop_is_EBH_reduction;
    LGMdouble EBH_reduction_parameter;
    bool growthloop_is_radiation_use_efficiency;
    LGMdouble radiation_use_efficiency_parameter;
    LGMdouble ebh_final_value;   ///< If -EBHREDUCTION is set, ebh parameters reach this value
    int growthloop_ebh_mode;     ///< Which variable (Qin, Qabs, or rue*Qabs) runs EBH growth distribution
    bool growthloop_is_heightFun; ///< If length of stem apical shoot is derived from relative crown length
    LGMdouble Db_current, Db_previous, dDb;  ///< For change in base diameter, for -heightFun
  };
}//End namespace LignumForest
#endif
#include <GrowthLoopI.h>
