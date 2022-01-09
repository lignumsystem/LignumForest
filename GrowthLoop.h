#ifndef GROWTH_LOOP_H
#define GROWTH_LOOP_H
#include <cmath>
#include <cstdio>
#include <fstream>
#include <utility>
#include <mathsym.h>
#include <Point.h>
#include <PositionVector.h>
#include <Bisection.h>
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
///\class CollectP
///\brief Collect photosynthates from the tree
///Use with Accumulate algorithm in stl-lignum 
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


///\class GrowthLoop GrowthLoop.h
///\brief Implement the LignumForest growth loop
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
     pairwise_self(false), eero(false), g_fun_var(0.0), g_fun_varies(false),
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

  void usage()const;
  void checkCommandLine(int argc, char** argv)const;
  void parseCommandLine(int argc, char** argv);
  void resolveCommandLineAttributes();

  void setVerbose(bool wordy){verbose=wordy;}
  bool isVerbose()const{return verbose;}
  void initRan3(int init)const{
    int neg_init=-abs(init);
    ran3(&neg_init);
  }
  void createTrees();
  void initializeTrees();
  void initializeVoxelSpace();
  void initializeFunctions();
  void initializeGrowthLoop();



  void increaseXi(int& year);
  void setTreeLocations();
  void photosynthesis(TREE& t);
  void respiration(TREE& t);
  void collectDataBeforeGrowth(TREE& t,unsigned int i);
  void treeAging(TREE& t);
  double collectSapwoodMass(TREE& t);
  void setSapwoodDemandAtJunction(TREE& t);
  bool allocation(TREE& t,bool verbose);
  void writeOutput(TREE& t,unsigned int tree_n,int iteration);
  void writeSensitivityAnalysisData(TREE& t);
  void writeCrownLimitData(TREE& t,int iteration);
  void writeVoxels(TREE& t);
  void writeTreeToXMLFile(TREE& t,int age,int interval)const;
  void writeBranchInformation(TREE& t,const string& file)const;
  void writeProductionBalance(TREE& t,const string& file)const;
  //vertical distribution of fip
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
  void setVoxelSpaceAndBorderForest();
  void calculateRadiation();
  StandDescriptor<TREE>& getStand() {return stand;}
  vector<TREE*>& getTrees() {return  vtree;}
  void increaseXi();
  void photosynthesisAndRespiration();
  void evaluateStandVariables();
  void createNewSegments();
  void allocationAndGrowth();
  void output();
  int getIterations() {return iterations;}
  int getNumberOfTrees() {return no_trees;}
  void setYear(const int& y) {year = y;}
  void setHPrev(){
    for (unsigned int i = 0; i < (unsigned int)no_trees; i++)
      h_prev[(int)i] = GetValue(*vtree[i],LGAH);
  }
private:
  vector<TREE*> vtree;//vector of trees
  vector<LSYSTEM*> vlsystem;//vector of L-systems, one for each tree
  vector<pair<double,double> > locations;
  vector<ofstream*> vdatafile;
  VoxelSpace *vs;//The voxel space
  StandDescriptor<TREE> stand;
  BorderForest border_forest;
  bool verbose;
  bool bracket_verbose;
  int iterations;
  int start_voxel_calculation;
  int num_parts;
  int interval;//write interval in output
  int init_density;//Initial density in self thinning scheme
  int nsegment;//Number of segments
  double tree_distance;//Minimum distance between trees in meters
  int  hw_start;//Starting year of heartwood build up
  double light_method;//Mandatory method for radiation assessment
  vector<double> wsapwood; //Sapwood mass  before new growth, BEFORE
  //sapwood senescence
  vector<double> wfoliage; //Foliage  mass before new  growth
  vector<double> wroot; //Root mass before new growth
  vector<double> ws_after_senescence; //Sapwood   mass   before  new
  //growth, AFTER sapwood senescence
  double ws2; //Old segment sapwood mass
  double ws3; //Sapwood mass allocated to diameter growth (ws2-ws1)
  double ws4; //Sapwood mass in the new segments
  double wr1; //Root mass required by new foliage
  double wfnew; // New foliage mass
  double wfaftergrowth;//Foliage after new growth
  double wsaftergrowth;//Sapwood after new growth
  double qintop1; //Qin at the top of the tree
  double cvol; //crown volume
  double axisvol;//Main stem volume
  double qabs; //Qabs
  double treeAf; //Foliage area of the tree
  double ASeg0; //Surface area of new segments
  double wood; //Wood mass of tghe tree
  double wstem; //wood mass in the main axis
  double wbranches; //wood mass in branches
  double stem_sw; //sapwood mass in the main stem
  double treeP; //Photosynthetic production of the tree. WARNING: be
  //careful when  using this  variable.  It is  set in
  //allocationAndGrowth() (called by allocation()) for
  //each  tree.  After  the vector  of trees  has been
  //passed over by allocation(), treeP holds the value
  //of the last tree in the vector
  double treeM; //Respiration of the tree. WARNING: the same warning
  //as for treeP.
  double lambda; //lambda value after allocation (P-M-G(lambda) = 0)
  summing bs;//Mean branch lengths
  DCLData dcl;//Diameter and heigth at the crown base.
  CrownVolume<TS,BUD> cv;//crown_volume
  MainAxisVolume<TS,BUD> mav; //Main stem volume
  string metafile;
  string voxelfile;
  string xmlfile; //XML file where the tree can be saved and restored from
  string datafile;//Text file to save tree data each time step
  string fipfile;
  bool to_file;//Write outputfile
  bool sensitivity_analysis;//Data for sensitivity analysis
  bool crown_limit_data;
  bool writevoxels;
  bool increase_xi;//Increase LGPxi value
  int xi_start; //Starting year of Xi increase
  bool self_thinning;//self thinning harvesting scheme
  ParametricCurve K;//The 'K' function 
  ParametricCurve stems_ha;//density of the forest
  ParametricCurve fdensity_ha;//density of the border forest
  Sensitivity<TS,BUD> sensitivity; //sensitivity analysis file
  bool generate_locations;
  string location_file;
  ifstream location_stream;
  int no_trees;
  bool noWoodVoxel;
  bool wood_voxel;
  int year;
  unsigned int target_tree;   //one of the trees can be identified as target tree
  bool write_output;
  ofstream* stand_output;
  StandDescriptor<TREE> center_stand;
  ofstream* cstand_output;
  bool evaluate_border_forest;
  LGMdouble k_border_conifer;        //extinction coeffient in border forest for conifer foliage
  vector<int> no_h;
  vector<double> h_prev;
  LGMdouble p0_var;    //Random variation of p0 +- max <value> per cent from the value in Tree.txt
  string stand_file;
  string cstand_file;
  LGMdouble seg_len_var;
  //===============  24.10.2019
  bool pairwise_self;
  bool eero;
  double g_fun_var;
  bool g_fun_varies;
  bool random_branch_angle;
  double ba_variation;
  //=================== 30.9.2021
  bool growthloop_is_EBH;
  bool growthloop_is_EBH1;
  double growthloop_EBH1_value;
};

#endif
#include <GrowthLoopI.h>
