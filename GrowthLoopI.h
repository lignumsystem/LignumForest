#ifndef GROWTHLOOPI_H
#define GROWTHLOOPI_H
#include <string>
#include <glob.h>
#include <LignumForestGlobals.h>
///\file  GrowthLoopI.h
///
///GrowthLoopI.h runs completely in voxel space.
///Changes marked with **run-voxel**.

namespace Pine{

  extern int mode;  ///<Declared in L-system. 

}

using namespace Pine;

using namespace LignumForest;
// extern int LignumForest::ran3_seed; ///<Initialized in GrowthLoopI.h
// extern double LignumForest::H_0_ini; ///<Initialized in GrowthLoopI.h
// extern double LignumForest::H_var_ini;///<Initialized in GrowthLoopI.h
// extern int LignumForest::n_buds_ini_min;///<Initialized in GrowthLoopI.h
// extern int LignumForest::n_buds_ini_max;///<Initialized in GrowthLoopI.h
// extern double L_age; ///<Initialized in GrowthLoopI.h
// extern double LignumForest::rel_bud;///<Initialized in GrowthLoopI.h
// extern bool LignumForest::bud_variation;///<Initialized in GrowthLoopI.h
// extern double LignumForest::branch_angle;///<Initialized in GrowthLoopI.h
// extern ParametricCurve LignumForest::bud_view_f;///<Initialized in GrowthLoopI.h
// extern bool LignumForest::is_bud_view_function;///<Reinitialized in GrowthLoopI.h.
// extern double LignumForest::global_hcb;///<Initialized in GrowthLoopI.h
// extern double L_H;///<Initialized in GrowthLoopI.h
namespace LignumForest{                                 
  template <class TREE,class TS, class BUD, class LSYSTEM>
  int  CreateTreeXMLDataSet(const GrowthLoop<TREE,TS,BUD,LSYSTEM>& gl, LGMHDF5File& hdf5_file,const string& group_name, const int interval)
  {
    const vector<TREE*>& vt = gl.getTrees();
    TREE* t0 = vt[0];
    int age = static_cast<int>(GetValue(*t0,LGAage));
    if (age % interval == 0){
      XMLDomTreeWriter<TS,BUD> writer;
      string age_group = to_string(age)+"/";
      string new_group = group_name+age_group;
      //Groups must be created one by one
      hdf5_file.createGroup(new_group);
      string tree_xml;
      for (unsigned int i=0; i < vt.size(); i++){
	TREE* ti = vt[i];
	int tree_id = static_cast<int>(GetValue(*ti,TreeId));
	string id = to_string(tree_id);
	string dset_name = "Tree_"+id;
	tree_xml = writer.xmlToString(*ti);
	hdf5_file.createDataSet(new_group+dset_name,tree_xml);
	tree_xml.clear();
      }
    }
    return 0;
  }

  template<class TREE, class TS, class BUD, class LSYSTEM>
  GrowthLoop<TREE,TS,BUD,LSYSTEM>::~GrowthLoop()
  {
    for (unsigned int i = 0; i < vtree.size(); i++){
      delete vtree[i];
    }

    for (unsigned int i = 0; i < vlsystem.size(); i++){
      delete vlsystem[i];
    }

    for (unsigned int i = 0; i < vdatafile.size(); i++){
      vdatafile[i]->close();
      delete vdatafile[i];
    }

    if (stand_output){
      stand_output->close();
      delete stand_output;
    }
  }
  
  template<class TREE, class TS, class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::insertMetaFiles(const string& regexp)
  {
    if (verbose){
      cout << "INSERT META FILES (Glob expressions): " << regexp << endl;
    }
    glob_t glob_result;
    glob(regexp.c_str(),GLOB_TILDE|GLOB_BRACE,NULL,&glob_result);
    for (unsigned int i=0; i < glob_result.gl_pathc; ++i){
      string fname = glob_result.gl_pathv[i];
      metafile_q.push_back(fname);
    }
    sort(metafile_q.begin(),metafile_q.end(),less<string>());
    if (verbose){
      for (unsigned int i=0; i < metafile_q.size();i++){
	string fname = metafile_q[i]; 
	cout << "METAFILE: " << i << " " << fname <<endl;
      }
    }
  }
  
  template<class TREE, class TS, class BUD, class LSYSTEM>
  string GrowthLoop<TREE,TS,BUD,LSYSTEM>::popMetaFile()
  {
    string s = metafile_q.front();
    metafile_q.pop_front();
    return s;
  }

  template<class TREE, class TS, class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::insertModeChangeYears(const string& years)
  {
    stringstream year_ss(years);
    string year_string;
    while (getline(year_ss,year_string,',')){
      int year = atoi(year_string.c_str());
      modechangeyears_q.push_back(year);
    }
    //Sort to ascending order
    sort(modechangeyears_q.begin(),modechangeyears_q.end(),less<int>());
    if (verbose){
      cout << "Mode change years: " << flush;
      for (unsigned int i=0; i < modechangeyears_q.size();i++){
	int year = modechangeyears_q[i]; 
	cout << year <<" "<<flush;
      }
      cout <<endl;
    }
  }

   template<class TREE, class TS, class BUD, class LSYSTEM>
   int GrowthLoop<TREE,TS,BUD,LSYSTEM>::popModeChangeYear()
   {
     int year = modechangeyears_q.front();
     modechangeyears_q.pop_front();
     return year;
   }
  
  /// \snippet{lineno} GrowthLoopI.h Usagex
  // [Usagex]
  template<class TREE, class TS, class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::usage()const
  {
    cout << "Usage:  ./lig-forest -iter <value>  -metafile <file>  -hdf5 <file> -voxelspace <file>" <<endl;
    cout << "[-numParts <parts>]  [-treeDist <dist>] [-hw <hw_start>]" <<endl;
    cout << "[-toFile <filename> OBSOLETE] [-xml <filename> OBSOLETE] [-writeVoxels] [-sensitivity <filename>] " <<endl;
    cout << "[-fipdistrib <filename>] [-writeInterval interval]" << endl;
    cout << "[-seed <num>] [-increaseXi <value>] [-targetTree <num>] " <<endl;
    cout << "[-treeFile <filename> OBSOLETE] [-generateLocations  <num>] [-woodVoxel] [-treeLocations <file>]" << endl;
    cout << "[-writeOutput OBSOLETE] [-verbose] [-bracketVerbose] [-noBorderForest] [-seed <value>] [-kBorderConifer <value>]"  << endl;
    cout << "[-H_0_ini <value>] [-H_var_ini <value>] [-n_buds_ini_min <num>] [-n_buds_ini_max <vlaue>]" << endl;
    cout << "[-p0Var <value>] [-segLenVar <value>] [-pairwiseSelf] [-budVariation <value>] [-eero]" << endl;
    cout << "[-gFunVar <value>] [-branchAngleVar <value>]" << endl;
    cout << "[-space0] [-space1] [-space2] [-adHoc]" << endl;
    cout << "[-budViewFunction] [-EBH] -EBH1 <value>]" << endl;
    cout << "[-space2Distance <Value>] [-EBHREDUCTION <value>] [-EBHFINAL <value>] [-EBHInput <int>] [-RUE <value>]" << endl;
    cout << "[-modeChange <year1,year2,...,yearN] [-architecureChange <year>]" << endl;
    cout << "[-fsapwdown <file>]" << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "-iter Number of years to simulate" << endl;
    cout << "-metafile <globexpr> Glob expression (as for Unix command line) for MetaFiles." << endl 
	 << "                     Files (traditionally called Metafile*.txt) containg file" <<endl
	 << "                     locations for Tree parameters, Firmament configuration" << endl
	 << "                     and Tree functions." << endl;
    cout << "                     Glob expression can define alternative matchings with curly braces," << endl
         << "                     e.g. {MetaFile,MetaFile1}.txt matches both MetaFile.txt and MetaFile1.txt." << endl;
    cout << "                     The MetaFiles defined are consumed in alphabetical order" <<endl;
    cout << "-hdf5 HDF5 file for simulation results. Trees as XML strings are in the HDF5 file with TREEXML_PREFIX prefix. See -writeInterval" << endl;
    cout << "-writeInterval <number> Save trees in XML format in every `number` of years" << endl;   
    cout << "-generateLocations <num>  In this case <num> trees will be generated to random locations. If this" << endl;
    cout << "          is not on, tree locations will be read from file Treelocations.txt. This file can be changed" << endl;
    cout << "          by -treeLocations <file>. If location file is not found program stops." << endl;
    cout << "-woodVoxel                If woody parts are dumped to voxels (default = true)" << endl;
    cout << "-treeDist <dist>          Minimum distance between two trees (default = 0), works only with -generateLocations." << endl;  
    cout << "-numParts <parts>         Segments can be dumped to voxels by parts (i.e. they may belong to different voxels," << endl;
    cout << "                          default = 1)" << endl;
    cout << "-hw <hw_start>            Starting year of formation of heartwood Default = 15 years" << endl;
    cout << "-increaseXi <value>       If parameter ksi (fraction of heartwood in new segments) increases, <value> = starting year (default = 15)"
	 << endl;
    cout << "-xiIncrement <value>      The parameter xi increases by the rate of 0.1/<value>. Default is 0.1/25.0." <<endl;
    cout << "                          This means increase = 0.004/year after year 15 up to value 0.85 (as in FPB 35: 964-975 2008 article)"
	 << endl;
    cout << "-targetTree <num>         Any one of the trees can be identified as target tree (default = 0)" << endl;
    cout << "-writeOutput               Most of the things are written to their respctive file at -writeInterval interval (default false)" << endl;
    cout << "-verbose                  Output of progress of run if set (default = not set)." << endl;
    cout << "-bracketVerbose             If set, iteration information is printed out in allocation (default = not set)." << endl;
    cout << "-noBorderForest             No border forest around the stand (default = there is border forest)"
	 << endl;
    cout << "-seed <value>               seed for random number generator." << endl;
    cout << "-kBorderConifer <value>       Extinction coefficient for conifer foliage in border forest (default = 0.14)" << endl;
    cout << "-H_0_ini, -H_var_ini        For variation of initial heights (defaults = 0.3 and 0.0)" << endl;
    cout << "-n_buds_ini_min, -n_buds_ini_max  For variation of initial number of buds (defaults = 4 and 4)" << endl;
    cout << "-p0Var <value>                Random variation of p0 +- max <value> per cent from the value in Tree.txt" << endl;
    cout << "-segLenVar         Random variation of length of new segments around Lnew, per cent" << endl;
    cout << "-pairwiseSelf      Pairwise radiation calculation for the tree itself." << endl;
    cout << "-eero              For studying the relationships between a) leaf light climate and b) syncronized variation" << endl;
    cout << "                   in leaf nitrogen concentration, leaf mass per area and leaf longevity, shoot length and" << endl;
    cout << "                   shoot leaf area and c) axis thickness scaling from the base of the stem to the axis tips." << endl;
    cout << "EBH resource distn can be in use in two ways. Both are set by command line arguments." << endl;
    cout << "-EBH               means EBH is in use and values (of lambda parameter) are" << endl;
    cout << "                   specified for all Gravelius orders (function SPEBHF, ScotsPine.h, function" << endl;
    cout << "                   file is specified the constructor of the tree)." << endl;
    cout << "-EBH1 <value>      means that EBH is in use and one value <value> is used" << endl;
    cout << "                   for all Gravelius orders. Option -EBH1 <value> overrides option -EBH" << endl;
    cout << "                   (EBH is set as SPis_EBH (Scots Pine Parameter Double SPPD), thus 0 == false, 1 == true)" << endl;
    cout << "                   EBH is according to W. Palubicki and K. Horel and S. Longay and" << endl;
    cout << "                   A. Runions and B. Lane and R. Mech and P. Prusinkiewicz. 2009." << endl;
    cout << "                   Self-organizing tree models for image synthesis ACM Transactions on Graphics 28 58:1-10." << endl;
    cout << "-EBHREDUCTION <value> If values of EBH parameters for all orders are reduced or increased  as" << endl;
    cout << "                      new_value = <value> * prev_value in each year after year 20. The goal of change" << endl;
    cout << "                   can be set by -EBHFINAL. The default value is 0.5"<< endl;
    cout << "-EBHFINAL <value>  Sets the goal value of -EBHREDUCTION" << endl;
    cout << "-EBHInput <int>    Changes the variable that runs EBH. If int == 1, it is Qabs, if int == 2, it is rue*Qabs," << endl;
    cout << "                   any other value or missing -EBHInput means Qin runs EBH." << endl;
    cout << "-RUE <value>       The radiation use effeciency (rue) varies as a function of TreeSegments initial radiation" << endl;
    cout << "                   conditions. Photosynthetic production of TreeSegment = rue * LGApr * Qabs. <value> = degree of" << endl;
    cout << "                   increase of rue as a function of shadiness (0 < <value> < 2)." << endl;
    cout << "-modeChange <year1,year2,..,yearN> Comma separated list of years when to apply new MetaFile." <<endl;      
    cout << "-architectureChange <year>  Change the braching pattern in L-system after <year> in simulation." <<endl;
    cout << "-Lmaxturn          Turn angle in degrees in the side branches when architecture change is on (default 80 degrees)." <<endl;
    cout << "-fsapwdown <file>  Part of the sapwood going down in a tree as a function of Gravelius order." <<endl;
    cout << "-butt_swell_coeff <value>  Adjustment coefficient (between 0 and 1) for the butt swell model." <<endl;
    cout << "-butt_swell_start <value>  Tree age to start butt swell." <<endl;
    cout << "-terminate_buds    Terminate buds grown out of VoxelSpace (set status DEAD)" <<endl;
    cout << endl;
  }
  // [Usagex]

  template<class TREE, class TS, class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::checkCommandLine(int argc, char** argv)const
  {
    if (argc < 5){
      cout << "Five (5) mandatory command line arguments are required!" << endl << endl;
      usage();
      exit(0);
    }
    ///+ Mandatory argument -iter
    else if (CheckCommandLine(argc,argv,"-iter") == false){
      cout << "Mandatory -iter <num> option missing" << endl;
      exit(0);
    }
    ///+ Mandatory argument -metafile
    else if (CheckCommandLine(argc,argv,"-metafile") == false){
      cout << "Mandatory -metafile <regexp> option missing" << endl;
      exit(0);
    }
    ///+ Mandatory argument -voxelspace
    else if (CheckCommandLine(argc,argv,"-voxelspace") == false){
      cout << "Mandatory -voxelspace <VoxelSpace.txt> option missing" << endl;
      exit(0);
    }
    ///+ Mandatory argument -hdf5
    else if (CheckCommandLine(argc,argv,"-hdf5") == false){
      cout << "Mandatory -hdf5 <HDF5file.h5> option missing" <<endl;
      exit(0);
    }
    ///+ Mandatory argument -writeInterval
    else if (CheckCommandLine(argc,argv,"-writeInterval") == false){
      cout << "Mandatory -writeInterval <number> option missing" <<endl;
      exit(0);
    }
    ///+ Mandatory argument -fsapwdown
    else if (CheckCommandLine(argc,argv,"-fsapwdown") == false){
      cout << "Mandatory -fsapwdown <file> option missing" <<endl;
      exit(0);
    }
    else if (verbose){
      cout << "Command line O.K." <<endl;
    } 
  }

 
  template<class TREE, class TS, class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::parseCommandLine(int argc, char** argv)
  {
    //Check verbose
    verbose = false;
    if (CheckCommandLine(argc,argv,"-verbose")){
      verbose = true;
    }

    if (verbose){
      cout << "parseCommandLine begin" <<endl;
    }
    ///---
    ///\par Check the command line for mandatory arguments.
    checkCommandLine(argc,argv);
    ///---
    ///\par Parse command line arguments
    string clarg;
    ///\par Parse arguments for iterations, MetaFile and VoxelSpace
    ///
    ///+ `-iter`, number of iterations (i.e. years)
    if (ParseCommandLine(argc,argv,"-iter", clarg)){
      iterations = atoi(clarg.c_str());
    }
    ///+ `-metafile`, Regular expression for MetaFiles containing actual parameter files etc.
    ///\sa metafile metafile_q
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-metafile", clarg)){
      metafile = clarg;
      insertMetaFiles(metafile);
    }
     ///---
    ///\par Parse growth mode change year
    clarg.clear();
    ///+ -modeChange Growth mode change years
    ///\sa LignumForest::is_mode_change GrowthLoop::insertModeChangeYears()
    if (ParseCommandLine(argc,argv,"-modeChange",clarg)){
      LignumForest::is_mode_change=true;
      insertModeChangeYears(clarg);
      if (verbose){
	cout << "MetaFiles " << metafile_q.size() <<endl;
	cout << "Mode change years " << modechangeyears_q.size() <<endl;
      }
      //Consistency check: There must be one MetaFile more than Growth mode change years
      if (metafile_q.size() != modechangeyears_q.size()+1){
	cout << "Growth mode change consistency error" <<endl;
	cout << "The number of MetaFiles " << metafile_q.size()
	     << " must be the number of mode change years " <<  modechangeyears_q.size() << " plus 1" << endl;
	exit(0);
      }
    }
    ///+ `-voxelspace`, voxel space definition \sa voxelfile
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-voxelspace", clarg)){
      voxelfile = clarg;
      std::ifstream f(voxelfile);
      if (!f.good()){
	cout << "VoxelSpace file " << voxelfile  << " does not exists or not accessible" <<endl;
       	cout << "Exiting program" << endl;
	exit(0);
      }
    }
    ///+ `-pairwiseSelf` in voxel space radiation regime
    ///\sa pairwise_self
    pairwise_self = false;
    if (CheckCommandLine(argc,argv,"-pairwiseSelf")) {
      pairwise_self = true;
    }
    ///
    
    ///---
    ///\par Parse output files
    ///
    ///+ `-toFile`, output to data file. This is  the base name
    ///each tree will add its coordinates to make unique files.
    ///\sa datafile
    ///\deprecated Using HDF5 files instead
    clarg.clear();
    to_file = ParseCommandLine(argc,argv,"-toFile", clarg);
    if (to_file){
      datafile = clarg;
    }
    ///+ `-cstandFile`, default *cstand.dat*.
    ///\deprecated Using HDF5 files instead
    cstand_file = "cstand-values.dat"; 
    clarg.clear();
    ParseCommandLine(argc,argv,"-cstandFile", clarg);
    cstand_file = clarg;
    ///+ `-standFile`, default *stand-values.dat*.
    ///\deprecated Using HDF5 files instead
    stand_file = "stand-values.dat"; 
    clarg.clear();
    if(ParseCommandLine(argc,argv,"-standFile", clarg))
      stand_file = clarg;
    ///+ `-sensitivity`, sensitivity analysis file, no default value.
    ///\sa sensitivity_analysis sensitivity
    clarg.clear();
    sensitivity_analysis = ParseCommandLine(argc,argv,"-sensitivity", clarg);
    if (sensitivity_analysis){
      sensitivity.printHeader(clarg);
    }
    ///+ `-crownLimitData`, write  crown limit data.
    crown_limit_data =  CheckCommandLine(argc,argv,"-crowmLimitData");
    ///+ `-xml`, XML file where the tree can be saved and restored from.
    ///\deprecated Using HDF5 files instead
    clarg.clear();
    ParseCommandLine(argc,argv,"-xml", clarg);
    xmlfile = clarg;
    ///+ `-fipdistrib`, the vertical distribution of fip from segments.
    ///\deprecated Using HDF5 instead
    clarg.clear();
    ParseCommandLine(argc,argv,"-fipdistrib", clarg);
    fipfile = clarg;
    ///+ `-writeVoxels`, boolean write voxels, the file will be named as `VoxelSpace-x-y-age.txt`, default `false`.
    if (CheckCommandLine(argc,argv,"-writeVoxels")){
      writevoxels = true;
    }
    ///+ `-writeOutput`, boolean write output in the first place, default `false`.
    ///\deprecated Using HDF5 files
    write_output = false;
    if (CheckCommandLine(argc,argv,"-writeOutput"))
      write_output = true;

    ///---
    ///\par Parse arguments to control and conduct simulation
    ///
    ///+ `-numParts`, Number of segment parts used to assess Qabs.
    num_parts = 1;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-numParts", clarg))
      num_parts = atoi(clarg.c_str());
    ///+ `-noWoodVoxel`, boolean no woody parts into voxels, default `true`. 
    wood_voxel = true;
    if (CheckCommandLine(argc,argv,"-noWoodVoxel"))
      wood_voxel = false;
    ///+ `-writeInterval`, write interval (timesteps) in output, default `1`.
    interval = 1;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-writeInterval", clarg))
      interval = atoi(clarg.c_str());
    ///+ `-increaseXi`, increase LGPxi from a start year. \sa increase_xi xi_start. 
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-increaseXi", clarg)){
      increase_xi = true;
      xi_start = atoi(clarg.c_str());
    }
    ///+ `-xiIncrement <value>, LGPxi increment as 0.1/<value>. Default 0.1/25.0.
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-xiIncrement",clarg)){
      xi_increment=atof(clarg.c_str());
    }
    ///+ `-hw`, the start of heartwood build up, default `15` years.
    clarg.clear();
    hw_start = 15;
    if (ParseCommandLine(argc,argv,"-hw", clarg))
      hw_start = atoi(clarg.c_str());
    ///+ `-bracketVerbose`, boolean set allocation of photosynthates verbose mode, default `false`.
    ///\sa bracket_verbose
    bracket_verbose = false;
    if (CheckCommandLine(argc,argv,"-bracketVerbose"))
      bracket_verbose = true;
    ///.
    
    ///---
    ///\par Parse arguments for forest generation.
    ///
    ///+ `-generateLocations`, if present the argument is the number of trees.
    ///\attention Either `-generateLocations` or `-treeLocations` must be present
    ///\sa no_trees, generate_locations
    generate_locations = false;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-generateLocations", clarg)) {
      generate_locations = true;
      no_trees = atoi(clarg.c_str());
    }
    else {
      clarg.clear();
      ///+ `-treeLocations`, use file for tree locations, default `Treelocations.txt`.
      ///\pre `generate_locations == false`;
      if (ParseCommandLine(argc,argv,"-treeLocations", clarg)) 
	location_file = clarg;
      else
	location_file = "Treelocations.txt";
    }
    ///+ `-treeDist`, set tree distance.
    ///\sa tree_distance.
    tree_distance = 0.0;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-treeDist", clarg))
      tree_distance = atof(clarg.c_str());
    ///+ `-targetTree`, set the target tree to be analysed.\sa target_tree
    ///\deprecated Using HDF5 files instead
    target_tree = 0;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-targetTree", clarg))
      target_tree = (unsigned int)atoi(clarg.c_str());
  
    ///+ `-noBordeForest`, use  homogenous border forest, default `true`.
    ///\sa evaluate_border_forest.
    evaluate_border_forest = true;
    if (CheckCommandLine(argc,argv,"-noBorderForest"))
      evaluate_border_forest = false;

    ///+ `-seed`, seed for ran3() random number generator, default `-123321`.
    ///\sa LignumForest::ran3_seed. 
    LignumForest::ran3_seed = -123231;
    int s_ini;
    if (ParseCommandLine(argc,argv,"-seed", clarg)){
      if (clarg.length() > 0){
	s_ini = atoi(clarg.c_str());
	LignumForest::ran3_seed = -abs(s_ini);
      }
    }
    ///+ `-kBorderConifer`, extinction coefficient for homogenous conifer border forest, default 0.14.
    ///\sa k_border_conifer.
    k_border_conifer = 0.14;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-kBorderConifer", clarg))
      k_border_conifer = atof(clarg.c_str());
    ///.
    
    ///---
    ///\par Parse parameters for experiments to adjust tree growth.
    ///
    ///+ `-H_0_ini`
    LignumForest::H_0_ini = 0.3;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-H_0_ini", clarg))
      LignumForest::H_0_ini = atof(clarg.c_str());
    ///+ `-H_var_ini`
    LignumForest::H_var_ini = 0.0;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-H_var_ini", clarg))
      LignumForest::H_var_ini = atof(clarg.c_str());
    ///+ `-n_buds_ini_min`
    LignumForest::n_buds_ini_min = 4;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-n_buds_ini_min", clarg))
      LignumForest::n_buds_ini_min  = atoi(clarg.c_str());
    ///+ `-n_buds_ini_max`
    ///\sa H_0_ini, H_var_ini, n_buds_ini_min, n_buds_ini_max
    LignumForest::n_buds_ini_max = 4;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-n_buds_ini_max", clarg))
      LignumForest::n_buds_ini_max = atoi(clarg.c_str());
    ///+ `-p0Var` Random variation photosynthetic effieciency
    ///\sa p0_var
    p0_var = 0.0;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-p0Var", clarg))
      p0_var = atof(clarg.c_str());
    ///+ `-segLenVar` Activate variation in segment length in percentage
    ///\sa The Lignum::LGPlen_random parameter for tree set in Tree parameter file
    ///\sa LignumForest::is_random_length
    clarg.clear();
    if (CheckCommandLine(argc,argv,"-segLenVar")){
      LignumForest::is_random_length=true;
    }
    ///+ `-budVariation` in pecenatage
    ///\sa LignumForest::rel_bud LignumForest::bud_variation
    LignumForest::rel_bud = 0.0;
    LignumForest::bud_variation = false;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-budVariation", clarg)) {
      LignumForest::bud_variation = true;
      LignumForest::rel_bud = atof(clarg.c_str())/100.0;
    }
    ///+ `-eero`, Eero Nikinmaa model
    ///\sa eero
    eero = false;
    if (CheckCommandLine(argc,argv,"-eero")) {
      eero = true;
    }
    ///+ `-gFunVar` \sa g_fun_varies, g_fun_var
    g_fun_var = 0.0;
    g_fun_varies = false;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-gFunVar", clarg)) {
      g_fun_var = atof(clarg.c_str())*PI_VALUE/100.0;
      g_fun_varies = true;
    }
    ///+ `-branchAngleVar`
    ///\sa random_branch_angle ba_variation
    random_branch_angle = false;
    if (ParseCommandLine(argc,argv,"-branchAngleVar", clarg)) {
      random_branch_angle = true;
      ba_variation = atof(clarg.c_str())*PI_VALUE/100.0;
    }
    ///+ `-adHoc`
    ///\sa LignumForest::is_adhoc
    LignumForest::is_adhoc = false;
    if(CheckCommandLine(argc,argv,"-adHoc")) {
      LignumForest::is_adhoc = true;
    }
    ///+ `-budViewFunction`
    ///\sa LignumForest::is_bud_view_function
    LignumForest::is_bud_view_function = false;
    if(CheckCommandLine(argc,argv,"-budViewFunction")) {
      LignumForest::is_bud_view_function = true;
    }
    ///.
    LignumForest::space0 = false;
    LignumForest::space1 = false;
    LignumForest::space2 = false;
    
    ///---
    ///\par Parse space occupancy modelling.
    ///
    ///+ `-space0`
    if(CheckCommandLine(argc,argv,"-space0")) {
      LignumForest::space0 = true;
    }
    ///+ `-space1`
    if(CheckCommandLine(argc,argv,"-space1")) {
      LignumForest::space1 = true;
    }
    ///+ `-space2`
    if(CheckCommandLine(argc,argv,"-space2")) {
      LignumForest::space2 = true;
      clarg.clear();
      ///+ `-space2Distance`
      ///\sa LignumForest::space0  LignumForest::space1  LignumForest::space2  LignumForest::space2_distance
      if (ParseCommandLine(argc,argv,"-space2Distance",clarg)) {
	LignumForest::space2_distance = atof(clarg.c_str());
      }
    }
    ///.

    ///---
    ///\par Parse EBH model.
    ///
    ///+ `-EBH`
    growthloop_is_EBH = false;
    if(CheckCommandLine(argc,argv,"-EBH")) {
      growthloop_is_EBH = true;
    }
    ///+ `-EBH1`
    ///\sa growthloop_is_EBH growth_loop_is_EBH1 growthloop_EBH1_value
    ///\note Parameters are set in GrowthLoop::initializeTrees()
    growthloop_is_EBH1 = false;
    growthloop_EBH1_value = 0.0;
    if (ParseCommandLine(argc,argv,"-EBH1",clarg)) {
      growthloop_is_EBH1 = true;
      growthloop_EBH1_value = atof(clarg.c_str());
      growthloop_is_EBH = true;
    }
    ///+ Write `ebh.fun` function file.
    ///\pre growthloop_is_EBH == *true* and growthloop_is_EBH1 == *true*
    ///\sa growthloop_is_EBH  growthloop_is_EBH1
    ///\attention File name `ebh.fun` is hard coded 
    ///\snippet{lineno} GrowthLoopI.h ebhfun
    ///\internal
    //[ebhfun]
    if(growthloop_is_EBH  && growthloop_is_EBH1) { //a bit nonelegant but may work, see initializeTrees()
      ofstream fout("ebh.fun", ofstream::trunc);
      fout << "#Extended Borchert-Honda lambda as a f. of Gravelius order" << endl;
      fout << "#Main axis is 1 first branch 2 etc." << endl;
      LGMdouble x = 1.0;
      for(int i = 0; i < 3; i++) {
	fout << x << " "  <<  growthloop_EBH1_value << endl;
	x += 1.0;
      }
      fout << "6.0 " <<  growthloop_EBH1_value << endl;
      fout.close();
    }
    //[ebhfun]
    ///\endinternal
    ///+ `-EBHREDUCTION` \sa growthloop_is_EBH_reduction EBH_reduction_parameter
    growthloop_is_EBH_reduction = false;
    EBH_reduction_parameter = 1.0;
    clarg.clear();
    if(ParseCommandLine(argc,argv,"-EBHREDUCTION", clarg)) {
      growthloop_is_EBH_reduction = true;
      EBH_reduction_parameter = atof(clarg.c_str());
    }
    ///+ `-EBHFINAL` \sa ebh_final_value
    ebh_final_value = 0.5;
    clarg.clear();
    if(ParseCommandLine(argc,argv,"-EBHFINAL", clarg)) {
      ebh_final_value = atof(clarg.c_str());
    }
    ///+ `-EBHInput` \sa LignumForest::ebh_mode
    LignumForest::ebh_mode = 0;
    clarg.clear();
    if(ParseCommandLine(argc,argv,"-EBHInput", clarg)) {
      LignumForest::ebh_mode  = atoi(clarg.c_str());
    }
    ///.

    ///---
    ///\par Parse radiation use efficiency
    ///
    ///+ -RUE, radiation use efficiency \sa growthloop_is_radiation_use_efficiency radiation_use_efficiency_parameter 
    growthloop_is_radiation_use_efficiency = false;
    radiation_use_efficiency_parameter = 0.0;
    clarg.clear();
    if(ParseCommandLine(argc,argv,"-RUE", clarg)) {
      growthloop_is_radiation_use_efficiency = true;
      radiation_use_efficiency_parameter = atof(clarg.c_str());
    }
    ///.
    
    ///---
    ///\par Parse architure change year
    clarg.clear();
    ///+ -architectureChange Architecture change year \sa Pine::is_architecture_change Pine::architecture_change_year
    if (ParseCommandLine(argc,argv,"-architectureChange",clarg)){
	Pine::is_architecture_change = true;
	Pine::architecture_change_year = atoi(clarg.c_str());
    }
    ///---
    ///\par Parse Sapwood down function
    clarg.clear();
    ///+ -fsapwdown Function file for determining pipe model to pass sapwood down from branches to the main stem
    ///and down to the base of the tree   
    if (ParseCommandLine(argc,argv,"-fsapwdown",clarg)){
      fsapwdownfile = clarg;
    }
    ///---
    ///\par Parse max branch turn angle after architecture change
    clarg.clear();
    ///+ -Lmaxturn Parse max turn angle after architecture change. The argument angle is in degrees. 
    if (ParseCommandLine(argc,argv,"-Lmaxturn",clarg)){
      //Change degrees to radians used
      double turn_angle = atof(clarg.c_str())*PI_VALUE/180.0;
      LignumForest::max_turn_in_architecture_change = turn_angle;
    }
    ///---
    ///\par Parse butt swell coefficent
    if (ParseCommandLine(argc,argv,"-butt_swell_coeff",clarg)){
      double bsc = atof(clarg.c_str());
      if (bsc < 0.0 || bsc > 1.0){
	cout << "Butt swell coefficient not in interval [0,1] " << bsc <<endl;
	cout << "Exiting program" << endl;
	exit(0);
      }
      LignumForest::butt_swell_coeff = bsc;
      cout << "Butt swell coefficient " <<  LignumForest::butt_swell_coeff << endl;
    }
    ///---
    ///\par Parse butt swell start year
    if (ParseCommandLine(argc,argv,"-butt_swell_start",clarg)){
      double bss = atoi(clarg.c_str());
      LignumForest::butt_swell_start = bss;
      cout << "Butt swell start year " << LignumForest::butt_swell_start <<endl;
    }
    ///---
    ///\par Check terminate escaped buds option
    if (CheckCommandLine(argc,argv,"-terminate_buds")){
      LignumForest::terminate_escaped_buds = true;
      cout << "Terminate escaped buds " << LignumForest::terminate_escaped_buds << endl;
    }

    ///---
    ///\par Limiting length of new shoots
    string iline;
    ifstream input_file(LignumForest::SEGMENT_LENGTH_LIMIT_FILE);
    if (!input_file.is_open()) {
      cout << "Could not open: " <<  LignumForest::SEGMENT_LENGTH_LIMIT_FILE << endl;
      exit(-1);
    }
    getline(input_file,iline);
    getline(input_file,iline);
    getline(input_file,iline);
    input_file >> LignumForest::length_limit_year;
    input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    input_file >> LignumForest::g1maxL >> LignumForest::g2maxL;
    input_file.close();

    if (verbose){
      cout << "parseCommandLine end" <<endl;
      printVariables();
    }
    
  }//End parseCommandLine()




  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::resolveCommandLineAttributes()
  {
    if(eero) {
      p0_var = 0.0;
      LignumForest::bud_variation = false;
    }

  }



  template<class TREE, class TS, class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::printVariables()const
  {
    cout << "Printing variables governing the growth" <<endl;
    cout << "Verbose " << verbose << endl
	 << "Iterations " << iterations << endl
	 << "Start voxel calculations " << start_voxel_calculation << endl
	 << "Number of segment parts " << num_parts << endl
	 << "Write interval " << interval << endl
	 << "Self thinning density " << init_density << endl
      //       << "Gap radius " << gap_radius << endl
      //       << "Target tree location " << x_coord << " " << y_coord << endl
	 << "Minimum tree distance " << tree_distance << endl
	 << "Heartwood build-up " << hw_start << endl
	 << "Light method " << light_method << endl
	 << "Meta file " <<  metafile << endl
	 << "Voxel file " << voxelfile << endl
	 << "XML file " << xmlfile << endl
	 << "FIP file " << fipfile << endl
	 << "Write voxels " << writevoxels << endl
	 << "Increase LGPxi " << increase_xi << endl
      //       << "Bounding box " << bounding_box << endl
      //       << "Target voxel extinction " << target_voxel_extinction << endl
      //       << "Expand Gap " << expand_gap << endl
	 << "Self thinning scheme " << self_thinning <<endl;
  }


  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::createTrees()
  {
    for (int i = 0; i < no_trees; i++){
      pair<double,double> p = locations[i];
      LSYSTEM *l = new LSYSTEM();
      ///\par Default functions
      ///Hard coded default functions that should exist in the run directory
      ///\internal
      ///\snippet{lineno} GrowthLoopI.h DFun
      // [DFun]
      TREE* t = new TREE(Point(p.first,p.second,0.0),PositionVector(0,0,1),
			 "sf.fun","fapical.fun","fgo.fun",
			 fsapwdownfile,"faf.fun","fna.fun", "fwd.fun",
			 "flr.fun","ebh.fun","bvf.fun");
      // [DFun]
      ///\endinternal
      //Set the TreeId
      SetValue(*t,TreeId,static_cast<double>(i));
      if (verbose){
	cout << "Created a tree at: " << p.first << " " << p.second <<endl;
      }
      vtree.push_back(t);
      vlsystem.push_back(l);

      no_h.push_back(0);
      h_prev.push_back(0.0);

      wsapwood.push_back(0.0);
      wfoliage.push_back(0.0);
      wroot.push_back(0.0);
      ws_after_senescence.push_back(0.0);

      if(to_file) {
	//Create, open and add data file for the tree
	ostringstream data_file;
	unsigned int  n = datafile.find_last_of('.');
	if ( n == 0){
	  //No suffix by the user
	  data_file << datafile << "-" << p.first <<"-"<<p.second << ".txt";
	}
	else{
	  string prefix = datafile.substr(0,n);
	  data_file << prefix  << "-" << p.first <<"-"<<p.second << ".txt";
	}
	ofstream* f = new ofstream(data_file.str().c_str(),ios_base::trunc);
	if (verbose){
	  cout << "Created data file " << data_file.str() << " for the  tree" <<endl;
	}
	//write header here
	*f << "year H QinTop RelQinTop Dbase Dbh Dcb Hcb P M WfAfter BDiam Blen WrAfter Mr MAbove "
	  "WsAfter MswBefore MfBefore VCrown NSeg Qabs RelQabs WfBefore PSpecific VStem ASurNweSeg iWfElong "
	  "iWswElong iWswSecG iWr iWswGrowth WwoodAbov WStem WBranch WswStem AswBase AswBh AswCb lambda " 
	  "Af" << endl;
	vdatafile.push_back(f);
      }   // if(to_file) ...
    }
  }


  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::resizeTreeDataMatrix()
  {
    int iter = getIterations();
    int ntrees = getNumberOfTrees();
    hdf5_tree_data.resize(iter+1,ntrees,TREE_DATA_COLUMN_NAMES.size());
    hdf5_tree_data.init(std::nan(""));
    hdf5_stand_data.resize(iter+1,STAND_DATA_COLUMN_NAMES.size());
    hdf5_stand_data.init(std::nan(""));
    hdf5_center_stand_data.resize(iter+1,STAND_DATA_COLUMN_NAMES.size());
    hdf5_center_stand_data.init(std::nan(""));
    hdf5_dead_tree_data.resize(iter+1,ntrees,DEAD_TREE_DATA_COLUMN_NAMES.size());
    hdf5_dead_tree_data.init(std::nan(""));
    lambdav.resize(static_cast<unsigned int>(ntrees),std::nan(""));
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  TMatrix2D<double> GrowthLoop<TREE,TS,BUD,LSYSTEM>::getHDF5TreeParameterData()
  {
    TMatrix2D<double> tree_parameters(1,TREE_PARAMETER_NAMES.size(),0.0);
    map<string,double> tree_pvalues;
    if (!vtree.empty()){
      TREE* t = vtree[0];
      tree_pvalues["LGPaf"] = GetValue(*t,LGPaf);
      tree_pvalues["LGPapical"] = GetValue(*t,LGPapical);
      tree_pvalues["LGPar"] = GetValue(*t,LGPar);
      tree_pvalues["LGPlen_random"] = GetValue(*t,LGPlen_random);
      tree_pvalues["LGPLmin"] = GetValue(*t,LGPLmin);
      tree_pvalues["LGPlr"] = GetValue(*t,LGPlr);
      tree_pvalues["LGPmf"] = GetValue(*t,LGPmf);
      tree_pvalues["LGPmr"] = GetValue(*t,LGPmr);
      tree_pvalues["LGPms"] = GetValue(*t,LGPms);
      tree_pvalues["LGPna"] = GetValue(*t,LGPna);
      tree_pvalues["LGPnl"] = GetValue(*t,LGPnl);
      tree_pvalues["LGPpr"] = GetValue(*t,LGPpr);
      tree_pvalues["LGPq"] = GetValue(*t,LGPq);
      tree_pvalues["LGPrhoW"] = GetValue(*t,LGPrhoW);
      tree_pvalues["LGPsf"] = GetValue(*t,LGPsf);
      tree_pvalues["LGPsr"] = GetValue(*t,LGPsr);
      tree_pvalues["LGPss"] = GetValue(*t,LGPss);
      tree_pvalues["LGPxi"] = GetValue(*t,LGPxi);
      tree_pvalues["LGPzbrentEpsilon"] = GetValue(*t,LGPzbrentEpsilon);
      for (unsigned int i=0; i < TREE_PARAMETER_NAMES.size();i++){
	tree_parameters[0][i] = tree_pvalues[TREE_PARAMETER_NAMES[i]];
      }
    }
    return tree_parameters;
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  TMatrix2D<double> GrowthLoop<TREE,TS,BUD,LSYSTEM>::getHDF5TreeFunctionData(const LGMF fn_enum)
  {
    const ParametricCurve& fn = GetFunction(*vtree[0],fn_enum);
    vector<double> v = fn.getVector();
    TMatrix2D<double> fn_data = getLignumFnData(v);
    return fn_data;
  }  

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::initializeTrees()
  {
    LGMVERBOSE vrb=QUIET;
    if (verbose)
      vrb = VERBOSE;

    ///---
    ///\par The main mandatory initialization of trees
    ///For each tree in tree vector
    ///+ Initialize tree with Lignum::InitializeTree
    ///+ Set PineTree::SPHwStart year
    ///+ Set Lignum::TreeQinMax
    ///+ Set Lignum::LGPlen_random, randon variation in segment length
    ///+ Set PineTree::SPis_EBH  0.0 (off).
    string f = popMetaFile();
    InitializeTree<TS,BUD> init(f,vrb);
    for (unsigned int i=0; i < vtree.size(); i++){
      ///\internal
      ///\snippet{lineno} GrowthLoopI.h TINIT
      // [TINIT]
      SetValue(*vtree[i],SPHwStart,(double)hw_start);
      init.initialize(*vtree[i]);
      SetValue(*vtree[i],TreeQinMax,GetFirmament(*vtree[i]).diffuseBallSensor());
      //cout << "TREE SEG_LEN_VAR " << GetValue(*vtree[i],LGPlen_random) << endl;
      SetValue(*vtree[i], SPis_EBH, 0.0);
      // [TINIT]
      ///\endinternal
    }
    
    ///---
    ///\par Optional EBH model
    ///If EBH model (GrowthLoop::growthloop_is_EBH): <br>
    ///For each tree in tree vector
    ///+ Read EBH function file. Hard coded *ebh.fun* file required.
    ///+ Set PineTree::SPis_EBH  1.0 (on).
    ///+ Install PineTree::SPEBHF function
    if (growthloop_is_EBH){
      ///\internal
      ///\snippet{lineno} GrowthLoopI.h EBHmod
      // [EBHmod]
      ParametricCurve ebh("ebh.fun");
      for (unsigned int i=0; i < vtree.size(); i++){
	SetValue(*vtree[i], SPis_EBH, 1.0);
	SetFunction(*vtree[i],ebh,SPEBHF);
      }
      // [EBHmod]
      ///\sa GrowthLoop::growthloop_is_EBH
      ///\endinternal
    }
    ///---
    ///\par Set branch angle
    ///Setting branch angle for all trees. <br>
    ///Optional hard coded random variation in branch angle.
    ///\internal
    ///\snippet{lineno} GrowthLoopI.h BAVAR
    // [BAVAR]
    LignumForest::branch_angle = 45.0 * PI_VALUE / 180.0;
    for (unsigned int i=0; i < vtree.size(); i++){
      if(random_branch_angle) {
	LignumForest::branch_angle *= 1.0 + ((ran3(&LignumForest::ran3_seed)-0.5)/0.5) * ba_variation; 
      }
      vtree[i]->setBranchAngle(LignumForest::branch_angle);
    }
    // [BAVAR]
    ///\sa GrowthLoop::random_branch_angle
    ///\sa ScotsPineTree::setBranchAngle LignumForest::branch_angle
    ///\endinternal

    ///---
    ///\par Bud View Function to L-system
    ///The global bud view function LignumForest::bud_view_f  use is controlled by the global variable
    ///LignumForest::is_bud_view_function. If LignumForest::is_bud_view_function == false
    ///then the  LignumForest::bud_view_f is ignored.
    ///\note Bud view function file is specified in the constructor of the tree.
    ///\internal
    ///\snippet{lineno} GrowthLoopI.h BVFunc
    // [BVFunc]
    LignumForest::bud_view_f = GetFunction(*vtree[0], SPBVF);
    // [BVFunc]
    ///\endinternal
    
    ///---
    ///\par Optional random variation in the photosynthetic parameter
    ///Hard coded adjustment for photosynthetic efficiency 
    if (p0_var){
      for (unsigned int i=0; i < vtree.size(); i++){
	///\internal
	///\snippet{lineno} GrowthLoopI.h P0VAR
	// [P0VAR]
	LGMdouble p0 = GetValue(*vtree[i],LGPpr);
	p0 *= 1.0 + ((ran3(&LignumForest::ran3_seed) - 0.5)/0.5)*p0_var/100.0;
	SetValue(*vtree[i],LGPpr,p0);
	// [P0VAR]
	///\endinternal
      }
    }
    
    ///---
    ///\par Optional Eero Nikinmaa model
    ///Adjust LUE, SLA and foliage longetivity. Hard coded *eero.par* file required.
    if(eero) {
      for (unsigned int i=0; i < vtree.size(); i++){
	ifstream eero_file("eero.par");
	if (!eero_file) {
	  cout << "eero.par is missing!" << endl;
	  exit(0);
	}
 
	string line;
	getline(eero_file,line);//The header
	getline(eero_file,line); //The data line
	LGMdouble LUE_factor = 0.0;
	LGMdouble SLA_factor = 1.0;
	LGMdouble fol_age_factor = 1.0;
	istringstream iss(line);
	iss >> LUE_factor >> SLA_factor >> fol_age_factor;
	eero_file.close();

	///\par Light-use efficiency LUE
	///Hard coded LUE adjustment with random effect for photosynthetic parameter.
	///LUE adjustment and foliage respiration Lignum::LGPmf are also correlated.
	///\internal
	///\snippet{lineno} GrowthLoopI.h LUEVAR
	// [LUEVAR]
	LGMdouble LUE = GetValue(*vtree[i],LGPpr);
	LGMdouble LUE_change = 1.0 + ((ran3(&LignumForest::ran3_seed)-0.5)/0.5) * LUE_factor;
	if(LUE_change < 0.0)
	  LUE_change = 0.0;
	SetValue(*vtree[i],LGPpr,LUE_change*LUE);
	LGMdouble mf = GetValue(*vtree[i],LGPmf);
	SetValue(*vtree[i],LGPmf,LUE_change*mf);
	// [LUEVAR]
	///\endinternal
	
	///\par SLA adjustment
	///Hard coded random effect adjusted SLA function for trees. SLA adjustment assumes
	///that the function is defined at argument values 0.0, 1.0, and 1.1.
	///Find the first random number so that \f$\mathit{SLA_{change}} \geq 0.0 \f$.
	///\internal
	///\snippet{lineno} GrowthLoopI.h SLAVAR
	// [SLAVAR]
	LGMdouble SLA_change = 1.0 + ((ran3(&LignumForest::ran3_seed)-0.5)/0.5) * SLA_factor;
	while(SLA_change < 0.0) {
	  SLA_change = 1.0 + ((ran3(&LignumForest::ran3_seed)-0.5)/0.5) * SLA_factor;
	}
	ParametricCurve SLA_fun = GetFunction(*vtree[i], SPFSF);
	ostringstream to_new_SLA_fun;
	to_new_SLA_fun << 0.0 << " " << SLA_change*SLA_fun(0.0) << " " << 1.0 << " "
		       << SLA_change*SLA_fun(1.0) << " " << 1.1 << " " << SLA_change*SLA_fun(1.1);
	int dummy = 0;
	ParametricCurve new_SLA_fun(to_new_SLA_fun.str(), dummy);
	SetFunction(*vtree[i], new_SLA_fun, SPFSF);
	// [SLAVAR]
	///\endinternal
	
	///\par Foliage mortality
	///Hard coded random adjustment \f$\mathit{FMc}\f$ to foliage mortality.<br>
	///Foliage mortality function \f$FM_{fun}(\mathit{year})\f$ assumes that the function<br>
	///has been defined at values \f$ \mathit{year} \in \left\{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0\right\}\f$.<br>
	///Foliage mortality or foliage longevity changes are obtained by stretching the argument values.
	///\internal
	///\snippet{lineno} GrowthLoopI.h FMC
	// [FMC]
	LGMdouble FMc = 1.0 + ran3(&LignumForest::ran3_seed)*fol_age_factor;
	if(FMc < 0.0)
	  FMc = 0.0;
	ParametricCurve FM_fun = GetFunction(*vtree[i], LGMFM);      
	ostringstream to_new_FM_fun;
	to_new_FM_fun << 0.0 << " " << 1.0 << " "
		      << FMc*1.0 << " " << FM_fun(1.0) <<  " "
		      << FMc*2.0 << " " << FM_fun(2.0) <<  " "
		      << FMc*3.0 << " " << FM_fun(3.0) <<  " "
		      << FMc*4.0 << " " << FM_fun(4.0) <<  " "
		      << FMc*5.0 << " " << FM_fun(5.0) <<  " "
		      << FMc*6.0 << " " << FM_fun(6.0) <<  " "
		      << FMc*6.0 + 5.0 << " " << 0.0;
	dummy = 0;
	ParametricCurve new_FM_fun(to_new_FM_fun.str(), dummy);
	SetFunction(*vtree[i], new_FM_fun, LGMFM);
	// [FMC]
	///\endinternal
      }
    }   //eero

    ///---
    ///\par Optional variability in segment length
    ///Hard coded random effect on segment length as Gravelius order function.<br>
    ///Assumes that Gravelius order function has been defined for Gravelius orders 0-7,
    ///i.e. \f$ \mathit{f(go)}: \left\{ 0,1,\ldots,7 \right\} \mapsto [0,1] \f$
    ///The orders 0,1 are not changed.
    if(g_fun_varies) {
      for (unsigned int i=0; i < vtree.size(); i++){
	///\internal
	///\snippet{lineno} GrowthLoopI.h GOVAR
	// [GOVAR]
	ParametricCurve G_fun = GetFunction(*vtree[i],SPFGO);
	//Random effect
	LGMdouble G_c = 1.0 + ((ran3(&LignumForest::ran3_seed)-0.5)/0.5) * g_fun_var;
	//New function values
	double v2 = G_c*G_fun(2.0);double v3 = G_c*G_fun(3.0);double v4 = G_c*G_fun(4.0);
	double v5 = G_c*G_fun(5.0);double v6 = G_c*G_fun(6.0);double v7 = G_c*G_fun(7.0);
	//Function range check
	if(v2 > 1.0) v2 = 1.0;if(v3 > 1.0) v3 = 1.0;if(v4 > 1.0) v4 = 1.0;
	if(v5 > 1.0) v5 = 1.0;if(v6 > 1.0) v6 = 1.0;if(v7 > 1.0) v7 = 1.0;
	ostringstream to_new_G_fun;
	to_new_G_fun << 0.0 << " " << G_fun(0.0) << " "<< 1.0 << " " << G_fun(1.0) << " "
		     << 2.0 << " " << v2 << " " << 3.0 << " " << v3 << " " << 4.0 << " " << v4 << " "
		     << 5.0 << " " << v5 << " " << 6.0 << " " << v6 << " " << 7.0 << " " << v7;
	int dummy = 0;
	ParametricCurve new_G_fun(to_new_G_fun.str(), dummy);
	SetFunction(*vtree[i], new_G_fun, SPFGO);
	// [GOVAR]
	///\sa GrowthLoop::g_fun_varies GrowthLoop::g_fun_var 
	///\endinternal
      }
    }     //if(g_fun_varies)
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::initializeVoxelSpace()
  {
    ifstream vf(voxelfile.c_str());
    LGMdouble vx,vy,vz; vx = vy = vz = 0.0;
    LGMdouble s1,s2,s3,b1,b2; s1 = s2 = s3 = 0.0; b1=b2=0.0;
    vf >> vx >> vy >> vz >> s1 >> s2 >> s3 >> b1 >> b2;
    //vx,vy,vz define the voxel space  dimensions, s1, s2, s3 define the
    //voxel box dimensions, b1 and b2 the border forest width
    if (verbose){
      cout << "Voxel Space: " << vx << " " << vy << " " << vz << " " 
	   << s1 << " " << s2 << " " << s3 << " " << b1 << " " << b2 <<endl;
    }
    vs = new VoxelSpace(Point(0,0,0),Point(vx,vy,vz),
			s1,s2,s3,
			static_cast<int>(vx/s1),static_cast<int>(vy/s2),static_cast<int>(vz/s3),
			GetFirmament(*vtree[0]));
    //LignumForest::TerminateEscapedBuds stores VoxelSpace and CenterStand dimensions.
    //Set the height (z-coordinate) of upper right corner to infinity, i.e. max supported value by compiler.
    //A bud is then checked against (x,y) coordinates 
    terminate_buds.resize(Point(0,0,0),Point(vx,vy,std::numeric_limits<double>::max()),b1,b2);
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::initializeEscapedBuds(int width)
  {
    //Origo
    Point p0 = vs->getLowerLeftCorner();
    //Diaginally opposite point, clcokwise upper third point
    Point p6= vs->getUpperRightCorner();
    //Clockwise bottom third point
    Point p2(p6.getX(),p6.getY(),p0.getZ());
    //Clockwise bottom fourth point
    Point p3(p6.getX(),p0.getY(),p0.getZ());
    //Width (x direction) of the voxel space
    double p0p3 = p0||p3;
    //Length (y direction) of the voxels space
    double p2p3 = p2||p3;
    //LignumForest::terminateEscapedBuds stores VoxelSpace and CenterStand dimensions.
    //Set the height (z-coordinate) of upper right corner to infinity, i.e. max supported value by compiler.
    //A bud is then checked against (x,y) coordinates 
    terminate_buds.resize(Point(0,0,0),Point(p0p3,p2p3,std::numeric_limits<double>::max()),width,width);
  }
  
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::initializeFunctions()
  {
    ifstream fk("K.fun");
    if (fk.good()){
      K.install("K.fun");
    }
    ifstream fstem("stemsha.fun");
    if (fstem.good()){
      stems_ha.install("stemsha.fun");
    }
    ifstream fdensity("fdensity.fun");
    if (fdensity.good()){
      fdensity_ha.install("fdensity.fun");
    }
    if (verbose){
      cout << "K() O.K: " << K.ok() << " stems_ha() O.K: " << stems_ha.ok() 
	   << " fdensity() O.K: " << fdensity_ha.ok() << endl;
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::setTreeLocations()
  {  
    if(generate_locations) {
      ifstream vf(voxelfile.c_str());
      LGMdouble vx,vy,vz; vx = vy = vz = 0;
      LGMdouble s1,s2,s3,b1,b2; s1 = s2 = s3 = 0.0; b1=b2=0.0;
      vf >> vx >> vy >> vz >> s1 >> s2 >> s3 >> b1 >> b2;
    
      ///\par Forest stand set-up
      ///The opposite corners *l* and *r* of the *GrowthLoop::stand* from the voxel space.
      ///\internal
      ///\snippet{lineno} GrowthLoopI.h FArea
      // [FArea]
      Point l(0.0, 0.0, 0.0);
      Point r(vx, vy, 0.0);
      stand.setLlCorner(l);
      stand.setUrCorner(r);
      stand.evaluateArea();
      // [FArea]
      ///\endinternal
      
      ///\par Center stand set-up
      ///Set up GrowthLoop::center_stand that is used to evaluate stand values
      ///without border effect. The corners *cl* and *cr* from the
      ///voxel space
      ///\internal
      ///\snippet{lineno} GrowthLoopI.h CArea
      // [CArea]
      Point cl(b1,b2,0.0);
      Point cr(vx-b1, vy-b2, 0.0);
      center_stand.setLlCorner(cl);
      center_stand.setUrCorner(cr);
      center_stand.evaluateArea();
      // [CArea]
      ///\endinternal
      
      ///\par Border Forest set-up
      ///The opposite corners for GrowthLoop::border_forest as in GrowthLoop::stand.
      ///The border forest is the assumed to extend to infinity around the
      ///forest stand. The height of the homogenous layer of foliage of the border forest
      ///is dynamic based on the uppermost branches of the trees.
      ///\internal
      ///\snippet{lineno} GrowthLoopI.h BArea
      // [BArea]
      border_forest.setCornerL(l);
      border_forest.setCornerR(r);
      // [BArea]
      ///\endinternal
      
      ///\par Generate trees in the random locations
      ///The GrowthLoop::no_trees is limited by the GrowthLoop::tree_distance (so called hard core)
      ///\note The area of the *gap* is 0. The use of *gap* provides another hard core. This is also
      ///for historic reasons, consistent with LignumForest::GenerateLocations use in Lig-Crobas project.
      ///\internal
      ///\snippet{lineno} GrowthLoopI.h GenLoc
      // [GenLoc]
      //ForestGap gap((x_coord,y_coord),gap_radius);
      ForestGap gap(pair<double,double>(0.0,0.0),0.0);
      int no_trees_0 = no_trees;
      GenerateLocations(no_trees,0.0,0.0,vx,vy,tree_distance,gap,locations);
      // [GenLoc]
      ///\endinternal
      //Insert the target tree location
      //    locations.insert(locations.begin(), pair<double,double>(x_coord,y_coord));
      if (verbose){
	cout << "Number of trees" << locations.size() <<endl 
	     << " Density/ha wanted: " << (double)no_trees_0/(vx*vy/10000.0)
	     << " Density/ha created: " << (double)no_trees/(vx*vy/10000.0) <<endl;
	cout << " Minimum tree distance: " << tree_distance <<endl; 
      }
    } //  if(generate_locations  ...)
    else {
      ifstream location_stream(location_file.c_str());
      if(!location_stream) {
	cout << "Could not open tree location file " << location_file << endl;
	exit(0);
      }

      string line;
      getline(location_stream,line);
      LGMdouble cx, cy, bx, by;
      location_stream >> cx >> cy >> bx >> by;
      getline(location_stream,line);
      Point l(0.0, 0.0, 0.0);
      Point r(cx, cy, 0.0);
      stand.setLlCorner(l);
      stand.setUrCorner(r);
      stand.evaluateArea();

      //Set the same for center_stand that is used to evaluate stand values
      //without border effect
      Point cl(bx,by,0.0);
      Point cr(cx-bx, cy-by, 0.0);
      center_stand.setLlCorner(cl);
      center_stand.setUrCorner(cr);
      center_stand.evaluateArea();

      border_forest.setCornerL(l);
      border_forest.setCornerR(r);

      getline(location_stream,line);
      no_trees = 0;
      bool stop = false;
      ///\par Generate trees from the file
      ///Set-up forest, center stand and border forest as with random tree locations.
      ///Read tree location from the file.
      ///\internal
      ///\snippet{lineno} GrowthLoopI.h TFile
      // [TFile]
      while(!stop) {
	LGMdouble x, y;
	location_stream >> x >> y;
	getline(location_stream,line);
	if(!location_stream.eof()) {
	  no_trees++;
	  locations.insert(locations.end(), pair<double,double>(x,y));
	}
	else
          stop = true;
      }
      // [TFile]
      ///\endinternal
      if(no_trees < 1) {
	cout << "Reading tree locations from " << location_file << " did not succeed." << endl;
	exit(0);
      }
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::initializeGrowthLoop()
  {
    for (unsigned int i = 0; i < (unsigned int)no_trees; i++){
      LSYSTEM* l = vlsystem[i];
      TREE* t = vtree[i]; 

      l->start();
      l->lstringToLignum(*t,1,PBDATA);

      double wf = 0.0;
      wf = Accumulate(*t,wf,CollectNewFoliageMass<TS,BUD>());
      //Initial root mass
      SetValue(*t,TreeWr,GetValue(*t,LGPar)*wf);
      //Calculate  the  LGAsf  for   newly  created  segments,  sf  in  P
      //Kaitaniemi data depens on segment length
      ForEach(*t,SetScotsPineSegmentSf());
    }
    if (verbose)
      cout << "Initialized " << vtree.size() << " trees" <<endl;

    ///Obselete part, using HDF5 files.
    //Find suitable place for this: (maybe all output related things together)
    stand_output = new ofstream(stand_file.c_str(), ios_base::trunc);
    *stand_output << "year density DbaseAv DbaseMi DbaseMa DbhAv DbhMi DbhMa HAv HMi HMa "
      "G Gcb Vstem LAI Wf meanHcb"
		  << endl;
    cstand_output = new ofstream(cstand_file.c_str(),ios_base::trunc);
    *cstand_output << "year density DbaseAv DbaseMi DbaseMa DbhAv DbhMi DbhMa HAv HMi HMa "
      "G Gcb Vstem LAI Wf meanHcb"
		   << endl;

  }


  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::increaseXi(int& year)
  {
    if (increase_xi && year >= xi_start){
      bool first = true;
      for (unsigned int i = 0; i < (unsigned int)no_trees; i++) {
	TREE* t = vtree[i];
	///\par Increase the value of LGPxi
	///Increase the share of heartwood \f$\xi\f$ in new segments in \p year:
	///\f[
	///\xi = \left\{
	///  \begin{array}{ll}
	///   \xi+0.1/25.0 &: \xi \leq 0.85 \\
	///   0.85 &: \mathit{otherwise}
	///  \end{array}
	///  \right.
	///\f]
	///The share of heartwood increases as the tree matures, grows older
	///until LGPxi has reached its maximum value.
	///\internal
	///\snippet{lineno} GrowthLoopI.h XI
	//[XI]
	double xii = GetValue(*t, LGPxi);
	double xii_previous = xii;
	//LGPXi Increment
	xii += 0.1/xi_increment;
	//Check for max value
	if(xii > MAX_XII) xii=MAX_XII;
	SetValue(*t,LGPxi, xii);
	//[XI]
	///\endinternal
	if (verbose && first) {
	  cout << "Increasing LGPxi from " << xii_previous << " to " << xii <<endl;
	  first = false;
	}
      }
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::growthModeChange(int year)
  {
    if (LignumForest::is_mode_change && (!modechangeyears_q.empty()) && (year == nextModeChangeYear())){
      string f = popMetaFile();
      popModeChangeYear();
      LGMVERBOSE verb = verbose ? VERBOSE : QUIET;
      InitializeTree<TS,BUD> init(f,verb);
      for (unsigned int i=0; i < vtree.size(); i++){
	init.initialize(*vtree[i]);
      }
    }
  }
  
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::photosynthesis(TREE& t)
  {
    //TreeCompartments (i.e. segments) do their photosynthesis
    ForEach(t, TreePhotosynthesis<TS,BUD>());
    //and then sum  the photosynthetic rates of all  segments and update
    //tree's P
    LGMdouble initPh = 0.0;
    LGMdouble sumPh =.0;
    sumPh = Accumulate(t, initPh, SumTreePhotosynthesis<TS,BUD>());
    SetValue(t, TreeP, sumPh);
  }


  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::respiration(TREE& t)
  {
    //TreeCompartments (i.e. segments)  do photosynthesis
    ForEach(t, TreeRespiration<TS,BUD>());
    //and then sum respiration rates of all segments and update tree M
    LGMdouble sumM = 0.0;
    sumM = Accumulate(t, sumM, SumTreeRespiration<TS,BUD>());
    //Then the root respiration
    sumM += GetValue(t,LGPmr)*GetValue(t,TreeWr);
    //The respiration is the sum of the above and below ground respiration
    SetValue(t, TreeM, sumM);
  }

  
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::collectDataBeforeGrowth(TREE& t,unsigned int i)
  {
    wroot[i] = GetValue(t,TreeWr);  //root mass
    //Collect foliage 
    double wf = 0.0;
    wf = Accumulate(t,wf,CollectFoliageMass<TS,BUD>());
    wfoliage[i] = wf;
    //Collect sapwood mass
    double ws = 0.0;
    ws = Accumulate(t,ws,CollectSapwoodMass<TS,BUD>());
    wsapwood[i] = ws;
  }

  template<class TREE, class TS, class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::collectDataAfterGrowth(const int year,bool collect_stand)
  {
    for (unsigned int tree_num = 0; tree_num < vtree.size(); tree_num++){
      map<string,double> tdafter;//collect descriptive data or metrics first into this dictionary 
      summing bs; // Branch summaries for mean branch.
      DCLData dcl; //Crown base: Diameter and height.
      TREE& t = *vtree[tree_num];
      bs = Accumulate(t,bs,Branchmeans());
      dcl = AccumulateDown(t,dcl,AddBranchWf(),DiameterCrownBase<TS,BUD>());
      double crown_volume = cv(t);//Crown volume
      Point tpoint = GetPoint(t);
      double tree_id =  GetValue(t,TreeId);
      tdafter["TreeId"] = tree_id;
      tdafter["X"] = tpoint.getX();
      tdafter["Y"] = tpoint.getY();
      tdafter["Z"] = tpoint.getZ();
      int nsegment=0;
      tdafter["TreeNseg"]=static_cast<double>(AccumulateDown(t,nsegment,CountTreeSegments<TS,BUD>()));
      tdafter["TreeCrownVol"] = cv(t);
      tdafter["TreeH"] = GetValue(t,LGAH);
      Pine::L_H = tdafter["TreeH"]; //global L_H
      tdafter["TreeDBase"] = GetValue(t,LGADbase);
      tdafter["TreeDbh"] = GetValue(t,LGADbh);
      tdafter["TreeDCrownBase"] = dcl.DCrownBase();
      tdafter["TreeHCrownBase"] = dcl.HCrownBase();
      LignumForest::global_hcb = dcl.HCrownBase(); //Global hcb
      tdafter["TreeAsBase"] = GetValue(t,LGAAsbase);
      tdafter["TreeAsDbh"] = GetValue(t,LGAAsDbh);
      tdafter["TreeAsCrownBase"] = dcl.ASwCrownBase();
      double tree_af = 0.0;
      tdafter["TreeAf"] = Accumulate(t,tree_af,CollectFoliageArea<TS,BUD>());
      tdafter["AxisVol"] = mav(t);
      tdafter["TreeP"] = GetValue(t,TreeP);
      tdafter["TreeM"] = GetValue(t,TreeM);
      ///`wroot` collected *before* new growth.   
      ///`wsapwood` and `wfoliage` collected *before* new growth.
      ///\snippet{lineno} GrowthLoopI.h Prev 
      ///\internal
      // [Prev]
      tdafter["Mr_prev"] = GetValue(t,LGPmr)*wroot[tree_num];
      tdafter["M_above"] = GetValue(t,TreeP)- GetValue(t,LGPmr)*wroot[tree_num];
      tdafter["Ms"] = GetValue(t,LGPms)*wsapwood[tree_num];
      tdafter["Mf"] = GetValue(t,LGPmf)*wfoliage[tree_num];
      // [Prev]
      ///\endinternal
      double wf=0.0;
      tdafter["Wf"] = Accumulate(t,wf,CollectFoliageMass<TS,BUD>());
      double wf_new=0.0;
      tdafter["Wf_new"] = Accumulate(t,wf_new,CollectNewFoliageMass<TS,BUD>());
      double ws_old=0.0;
      double ws=0.0;
      tdafter["Ws"] = Accumulate(t,ws,CollectSapwoodMass<TS,BUD>());
      tdafter["Ws_old"] = Accumulate(t,ws_old,CollectOldSegmentSapwoodMass<TS,BUD>());
      ///`ws_after_senescence` vector collected *before* new growth
      ///\snippet{lineno} GrowthLoopI.h Ws
      ///\internal
      // [Ws]
      tdafter["Ws_D_growth"] = tdafter["Ws_old"]-ws_after_senescence[tree_num];
      // [Ws]
      ///\endinternal
      double ws_new=0.0;
      tdafter["Ws_new"] = Accumulate(t,ws_new,CollectNewSegmentSapwoodMass<TS,BUD>());
      tdafter["Ws_D_growth+Ws_new"] = tdafter["Ws_D_growth"]+tdafter["Ws_new"];
      tdafter["Wr"] = GetValue(t,TreeWr);
      tdafter["Wr_new"] = GetValue(t,LGPar)*tdafter["Wf_new"];
      list<TreeCompartment<TS,BUD>*>& ls = GetTreeCompartmentList(GetAxis(t));
      //GetTopQin located in Pine.h
      double qin_top=0.0;
      tdafter["QinTop"] = accumulate(ls.begin(),ls.end(),qin_top,GetTopQin<TS,BUD>());
      tdafter["QinMax"] = GetValue(t,TreeQinMax);
      tdafter["QinTop/QinMax"] =tdafter["QinTop"]/tdafter["QinMax"];
      double qabs=0.0;
      tdafter["Qabs"] = Accumulate(t,qabs,CollectQabs<TS,BUD>());
      tdafter["Qabs/DiffBallSensor"] = tdafter["Qabs"]/(GetFirmament(t).diffuseBallSensor()*tdafter["TreeAf"]);
      tdafter["Wf_P"] = wfoliage[tree_num];
      // No photosynthesis when collecting intial data
      if (!std::isnan(tdafter["TreeP"]/tdafter["Wf_P"]))
	tdafter["TreeP/Wf_P"] = tdafter["TreeP"]/tdafter["Wf_P"];
      else
	tdafter["TreeP/Wf_P"] = 0;
      double aseg0 = 0.0;
      // SurfaceAreaOfNewSegments already implemented as a concrete ScotsPineSegment type.
      tdafter["ASeg0"] = Accumulate(t,aseg0,SurfaceAreaOfNewSegments());
      double wood = 0.0;
      tdafter["W"] = Accumulate(t,wood,CollectWoodMass<TS,BUD>());
      tdafter["Wstem"] = GetValue(t,LGAWstem);
      tdafter["Wbranch"] = tdafter["W"]-tdafter["Wstem"];
      double ws_stem=0.0;
      tdafter["Ws_stem"] = Accumulate(t,ws_stem,CollectStemSapwoodMass<TS,BUD>());
      // Branch data
      // Max branch in each quadrant 
      vector<double> quadrantv(4,0.0);
      CrownGroundArea<TS,BUD> cga(t);
      vector<double>& res = Accumulate(t,quadrantv,cga);
      // Max branch from quadrants
      vector<double>::iterator max_branch = Lignum::max_elmnt<vector<double>::iterator>(res.begin(),res.end());
      tdafter["MaxBranch"] = *max_branch;
      tdafter["MeanBranch_SumD^2"] = bs.d2;
      tdafter["MeanBranch_SumL"] = bs.lsum;
      tdafter["MeanBranch_SumD^2*L"] = bs.d2l;
      // Branches special cases when collecting initial data:no branches in initial trees
      if (!std::isnan(bs.d2l/bs.d2))
	tdafter["MeanBranch_SumD^2*L/SumD^2"] = bs.d2l/bs.d2;
      else
	tdafter["MeanBranch_SumD^2*L/SumD^2"] = 0.0;
      tdafter["MeanBranch_Nbranch"] = bs.n_br;
      if (!std::isnan(bs.lsum/bs.n_br))
	tdafter["MeanBranch_SumL/Nbranch"] = bs.lsum/bs.n_br;
      else
	tdafter["MeanBranch_SumL/Nbranch"] = 0.0;
      tdafter["lambda"] = lambdav[tree_id];
      //Single row for for the tree, index the dictionary with the column names vector
      //The data will be in right order 
      for (unsigned int i=0; i < TREE_DATA_COLUMN_NAMES.size(); i++){
	hdf5_tree_data[year][tree_id][i]=tdafter[TREE_DATA_COLUMN_NAMES[i]];
      }
    }//End for loop
    // Collect data for 2D arrays from forest stand
    if (collect_stand){
      map<string,double> sdafter;//stand data
      map<string,double> csdafter;//center stand data
      sdafter["Year"] = year;
      csdafter["Year"] = year;
      sdafter["StandArea"] = stand.getArea();
      csdafter["StandArea"] = center_stand.getArea();
      sdafter["N_trees"] = stand.getNoTrees();
      csdafter["N_trees"] = center_stand.getNoTrees();
      sdafter["10000*N_trees/StandArea"] = 10000.0*stand.getNoTrees()/stand.getArea();
      csdafter["10000*N_trees/StandArea"] = 10000.0*center_stand.getNoTrees()/stand.getArea();
      sdafter["Dbase_mean"] = stand.getMeanDbase();
      csdafter["Dbase_mean"] = center_stand.getMeanDbase();
      sdafter["Dbase_min"] = stand.getMinDbase();
      csdafter["Dbase_min"] = center_stand.getMinDbase();
      sdafter["Dbase_max"] = stand.getMaxDbase();
      csdafter["Dbase_max"] = center_stand.getMaxDbase();
      sdafter["Dbh_mean"] = stand.getMeanDbh();
      csdafter["Dbh_mean"] = center_stand.getMeanDbh();
      sdafter["Dbh_min"] = stand.getMinDBH();
      csdafter["Dbh_min"] = center_stand.getMinDBH();
      sdafter["Dbh_max"] = stand.getMaxDBH();
      csdafter["Dbh_max"] = center_stand.getMaxDBH();
      sdafter["H_mean"] = stand.getMeanHeight();
      csdafter["H_mean"] = center_stand.getMeanHeight();
      sdafter["H_min"] = stand.getMinH();
      csdafter["H_min"] = center_stand.getMinH();
      sdafter["H_max"] = stand.getMaxH();
      csdafter["H_max"] = center_stand.getMaxH();
      sdafter["StandBasalArea"] = stand.getStandBasalArea();
      csdafter["StandBasalArea"] = center_stand.getStandBasalArea();
      sdafter["StandBasalAreaCrownBase"] = stand.getBasalAreaAtCrownBase();
      csdafter["StandBasalAreaCrownBase"] = center_stand.getBasalAreaAtCrownBase();
      sdafter["StandStemVol"] = stand.getStemVolume();
      csdafter["StandStemVol"] = center_stand.getStemVolume();
      sdafter["LAI"] = stand.getLAI();
      csdafter["LAI"] = center_stand.getLAI();
      sdafter["Stand_Wf"] = stand.getWfMass(); 
      csdafter["Stand_Wf"] = center_stand.getWfMass();
      sdafter["CrownLimit_mean"] = stand.getMeanCrownLimit();
      csdafter["CrownLimit_mean"] = center_stand.getMeanCrownLimit();
      for (unsigned int i = 0; i < STAND_DATA_COLUMN_NAMES.size(); i++){
	hdf5_stand_data[year][i] = sdafter[STAND_DATA_COLUMN_NAMES[i]];
	hdf5_center_stand_data[year][i] = csdafter[STAND_DATA_COLUMN_NAMES[i]];
      }
    }
    ///\remark Reinitialize the `lambdav` for the next growth cycle.
    unsigned int size = lambdav.size();
    lambdav.resize(size,std::nan(""));
  }
  
  template<class TREE, class TS, class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::collectDeadTreeDataAfterGrowth(const int year)
  {
    for (auto it=dead_trees.begin(); it != dead_trees.end(); advance(it,1)){
      //Map data name to its value
      map<string,double> dead_tree_data;
      TREE* t = vtree[*it];
      double tree_id =  GetValue(*t,TreeId);
      Point p = GetPoint(*t);
      double axis_vol = mav(*t);
      dead_tree_data["TreeId"]=tree_id;
      dead_tree_data["X"] = p.getX();
      dead_tree_data["Y"] = p.getY();
      dead_tree_data["Z"] = p.getZ();
      dead_tree_data["AxisVol"]=axis_vol;
      //Retrieve data and insert to 3d data table 
      for (unsigned int i=0; i < DEAD_TREE_DATA_COLUMN_NAMES.size(); i++){
	hdf5_dead_tree_data[year][tree_id][i]=dead_tree_data[DEAD_TREE_DATA_COLUMN_NAMES[i]];
      }
    }
  }

  template<class TREE, class TS, class BUD, class LSYSTEM>
  void GrowthLoop<TREE,TS,BUD,LSYSTEM>::removeDeadTreesAllOver()
  {
    //Now, dead trees are removed from everywhere
    auto It = vtree.begin();
    auto Is = vlsystem.begin();
    auto Il = locations.begin();
    //Vectors for data collected before new growth must be updated 
    auto Iws = wsapwood.begin();
    auto Iwf = wfoliage.begin();
    auto Iwr = wroot.begin();
    auto Iws_after = ws_after_senescence.begin();
    auto If =  vdatafile.begin();
    bool also_vdatafile = false;          //vdatafile vector may be empty
    if(vdatafile.size() > 0)
      also_vdatafile = true;
    auto In = no_h.begin();
    auto Ih = h_prev.begin();

    auto I = dead_trees.begin();
    unsigned int previous = 0;
    while(I != dead_trees.end())
      {
	unsigned int this_advance = *I - previous;
	advance(It, this_advance);
	advance(Is, this_advance);
	advance(Il, this_advance);
	if(also_vdatafile)
	  advance(If, this_advance);
	advance(In, this_advance);
	advance(Ih, this_advance);
	//Vectors for data before new growth
	advance(Iws,this_advance);
	advance(Iwf,this_advance);
	advance(Iwr,this_advance);
	advance(Iws_after,this_advance);
	
	delete *It;   //locations were not created by new
	delete *Is;
	if(also_vdatafile)
	  delete *If;


	cout << "Dead tree " << GetValue(**It,TreeId) << " at location " << Il->first << " " << Il->second << " was deleted in year "
	     << year << endl;

	vtree.erase(It);
	vlsystem.erase(Is);
	locations.erase(Il);
	//Vectors  for data  before new  growth must  also be  udated to
	//maintain the integrity of  the positions (indices) of trees to
	//these vectors
	wsapwood.erase(Iws);
	wfoliage.erase(Iwf);
	wroot.erase(Iwr);
	ws_after_senescence.erase(Iws_after);
	if(also_vdatafile)
	  vdatafile.erase(If);
	no_h.erase(In);
	h_prev.erase(Ih);

	previous = *I + 1;

	I++;
	no_trees--;
      }
  }
  
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::collectVoxelSpaceData(const int year, const int interval)
  {
    vsdata.insertData(*vs,year,interval);
  }
  
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::collectSapwoodAfterSenescence(TREE& t, unsigned int i)
  {
    ws = collectSapwoodMass(t);
    ws_after_senescence[i] = ws;
  }

  ///\remark TreeAging takes care of senescence in segments (above ground part) and root mortality.
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::treeAging(TREE& t)
  {
    
    ForEach(t,TreeAging<TS,BUD>()); 
    SetValue(t,TreeWr, 
	     GetValue(t,TreeWr)-GetValue(t,LGPsr)*GetValue(t,TreeWr));
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  double GrowthLoop<TREE, TS,BUD,LSYSTEM>::collectSapwoodMass(TREE& t)
  {
    LGMdouble ws = 0.0;
    ws = Accumulate(t,ws,CollectSapwoodMass<TS,BUD>());
    return ws;
  }


  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::setSapwoodDemandAtJunction(TREE& t)
  {
    double alku = 1.0;    //= Gravelius order of main axis
    PropagateUp(t,alku,SetSapwoodDemandAtJunction());
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  bool GrowthLoop<TREE, TS,BUD,LSYSTEM>::allocation(TREE& t, bool verbose)
  {
    
    try{
      /// \par The growth function G
      ///
      ///The growth function \f$G\f$ is hard coded in GrowthLoop::allocation.
      ///Using LignumForest::PartialSapwoodAreaDown where the sapwood area
      ///is passed down as such  between segments that are  in the same  axis.
      ///Only the segments  of higher gravelius order require  less sapwood  in a branching point.
      ///
      ///Bisection given initial values for \f$ \lambda \f$ (0 and 10) iterates its value
      ///trying to find \f$ P - M - G(\lambda) = 0 \f$ with given accuracy (0.01).
      ///\internal
      ///\snippet{lineno} GrowthLoopI.h GFUNCTION
      // [GFUNCTION]
      DiameterGrowthData data;
      LGMGrowthAllocator2<TS,BUD,LignumForest::SetScotsPineSegmentLength,
			  LignumForest::PartialSapwoodAreaDown,LignumForest::ScotsPineDiameterGrowth2,DiameterGrowthData>
	G(t,data,PartialSapwoodAreaDown(GetFunction(t,SPSD)));
      Bisection(0.0,100.0,G,0.01,verbose); 
      double tree_id = GetValue(t,TreeId);
      lambdav[tree_id] = G.getL();
      // [GFUNCTION]
      ///\sa LignumForest::PartialSapwoodAreaDown
      ///\endinternal
    }
    catch (TreeGrowthAllocatorException e){
      cout << "P < M " << e.getP() << " " << e.getM() << endl;
      return false;
    }
    catch(BisectionBracketException e){
      cout << "Could not bracket " << e.getFa() << " " << e.getFb() << " "  << e.getFbl()  <<endl;
      cout << e.getA() << " "  << e.getB() << " " << e.getBl() <<endl;
      return false;
    }
    catch(BisectionMaxIterationException e){
      cout << "Bisection Max iterations " << e.getMaxIter() << " exceeded "<< endl;
      cout << e.getFa() << " " << e.getFc() << " " << e.getFb() <<endl;
      return false;
    }
    //Allocation with Bisection is fine 
    return true;
  }


  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::allocationAndGrowth()
  {
    dead_trees.clear();
    //Create new buds by making derive with mode == 1 (global variable in L system)
    Pine::mode = 1;
    
    for (unsigned int k = 0; k < (unsigned int)no_trees; k++){
      cout << "Allocation loop with tree: " << k << " Number of trees " << no_trees <<endl; 
      TREE* t = vtree[k];
      LSYSTEM* l = vlsystem[k];

      /// \internal
      /// Initialize calculation of thickness growth induced by adding new shoots.
      /// \sa SetSapwoodDemandAtJunction
      /// \snippet{lineno} GrowthLoopI.h PipeModel
      // [PipeModel]
      double alku = 1.0;    //= Gravelius order of main axis
      PropagateUp(*t,alku,SetSapwoodDemandAtJunction());
      // [PipeModel]
      /// \remark This is for Pipe model calculations:
      /// \endinternal
      if(!allocation(*t,bracket_verbose)){
	//Iteration failed
	//Save position for deletion from lists for growing trees
	dead_trees.push_back(k);       //iteration failed, this tree is dead
	cout << "In GrowthLoop<TREE, TS,BUD,LSYSTEM>::allocationAndGrowth():" <<endl;
	cout << "   GrowthLoop<TREE, TS,BUD,LSYSTEM>:allocation failed()" <<endl;
	cout << "Pushed " << k << " to dead trees" <<endl;
	continue;
      }
      else if(GetValue(*t,LGAH) - h_prev[(int)k] < 0.001) {
	cout << L_age << " " << k << " " << no_h[(int)k] << endl;
	no_h[(int)k] += 1;
	if(no_h[(int)k] >= 3){
	  //Iteration failed
	  //Save position for deletion from lists for growing trees
	  dead_trees.push_back(k);
	  cout << "Tree height growth stagnant 3 successive years" <<endl;
	  cout << "Pushed " << k << " to dead trees" << endl;
	}
	continue;
      }

      //Calculate  LGAsf  for   newly  created  segments,  sf  in  P
      //Kaitaniemi data depens on segment length
      ForEach(*t,SetScotsPineSegmentSf());

      bool kill = false;
      PropagateUp(*t,kill,KillBudsAfterAllocation<TS,BUD>());

      //Now the lengths of the segments are such that G = P - M. Adjust the diameters
      // of old segments on the basis sapwood demand from above and those of new
      // segments on the basis of their length.
      DiameterGrowthData dgdata;
      AccumulateDown(*t,dgdata,PartialSapwoodAreaDown(GetFunction(*t,SPSD)),
		     ScotsPineDiameterGrowth2(LGMGROWTH));
      //Root growth
      double wfnew = 0.0;
      Accumulate(*t,wfnew,CollectNewFoliageMass<TS,BUD>());
      SetValue(*t,TreeWr, GetValue(*t,TreeWr) + GetValue(*t,LGPar)*wfnew);

      //Segment  dimensions   have  changed  in   allocation,  update
      //L-string
      l->lignumToLstring(*t,1,PBDATA);
      l->lstringToLignum(*t,1,PBDATA);

      //Pass the foliage mass of the mother segment to terminating buds and
      // synchronize with L-system
      double wftobuds = 0.0;
      PropagateUp(*t,wftobuds,ForwardWf<TS,BUD>());
      double lentobuds = 0.0;
      PropagateUp(*t,lentobuds,ForwardSegLen<TS,BUD>());
      l->lignumToLstring(*t,1,PBDATA);


      //============================================================================
      // Here is calculation of estimation of local needle area density for
      // using in estimation of bud fate
      //=============================================================================
      if (LignumForest::is_bud_view_function) {
	BoundingBox b1;
	bool foliage = true;
	FindCfBoundingBox<TS,BUD> fb1(foliage);
	b1 = Accumulate(*t, b1, fb1);

	Point ll = b1.getMin();
	Point ur = b1.getMax();

	VoxelSpace vs1(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
		       0.05,0.05,0.05,5,5,5,GetFirmament(*t));
	vs1.resize(ll, ur);
	vs1.reset();
	int num_parts = 5;
	DumpCfTree(vs1, *t, num_parts);

	LGMdouble cone_height = 0.5, cone_half_angle =  0.7; //= 40 degrees
	int no_points_on_rim = 12;
        
	SetBudViewFunctor sbvf(&vs1, cone_height, cone_half_angle, no_points_on_rim);
	ForEach(*t, sbvf);
      }
      LignumForest::branch_angle = t->getBranchAngle();        //this global variable goes to pine-em98.L
      l->derive();
      l->lstringToLignum(*t,1,PBDATA);
    }  //for (unsigned int k = 0; k < ...	  
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::output()
  {
    if(write_output) {
      *stand_output << year << " ";
      stand.writeOutput(*stand_output);
      *cstand_output << year << " ";
      center_stand.writeOutput(*cstand_output);
      //    for (unsigned int k = 0; k < (unsigned int)no_trees; k++){
      //      TREE* t = vtree[];
      TREE* t = vtree[target_tree];
      //Annual (time  step) output: tree, its position  in the vector
      //and iteration
      writeOutput(*t,target_tree,year);
      writeCrownLimitData(*t,year);
      //      writeVoxels(*t);
      writeTreeToXMLFile(*t,GetValue(*t,LGAage),interval);
      writeFip(*t,interval);
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::writeOutput(TREE& t, unsigned int tree_n,int iter)
  {
    if (to_file){
    
      //Calculate first values for output
      wfnew = 0.0;
      wfnew = Accumulate(t,wfnew,CollectNewFoliageMass<TS,BUD>());
      double treePcurrent = GetValue(t,TreeP);
      double treeMcurrent = GetValue(t,TreeM);
      //To plot how much is allocated to new growth (iWn), diameter growth (iWo)
      //and to roots (iWr) collect data
      ws2 = 0.0;
      ws2 = Accumulate(t,ws2,CollectOldSegmentSapwoodMass<TS,BUD>());
      ws3 = ws2-ws_after_senescence[tree_n];//diameter growth
      //ws4 together with wfnew is new growth in new segments above ground
      ws4 = 0.0;
      ws4 = Accumulate(t,ws4,CollectNewSegmentSapwoodMass<TS,BUD>());
      //This is needed for roots
      wr1 = GetValue(t,LGPar)*wfnew;

      Axis<TS,BUD>& axis = GetAxis(t);

      //Collect foliage after growth
      wfaftergrowth = 0.0;
      Accumulate(t,wfaftergrowth,CollectFoliageMass<TS,BUD>());
      //Collect sapwood mass after growth
      wsaftergrowth = 0.0;
      Accumulate(t,wsaftergrowth,CollectSapwoodMass<TS,BUD>());

      list<TreeCompartment<TS,BUD>*>& ls = GetTreeCompartmentList(axis);
      //Mean branch length
      bs.d2 = bs.d2l = bs.lsum = 0.0;
      bs.n_br = 0;
      bs = Accumulate(t, bs, Branchmeans() );

      //The Qin at the top
      qintop1 = 0.0;
      GetTopQin<TS,BUD> getTopQin1;
      qintop1 = accumulate(ls.begin(),ls.end(),qintop1, getTopQin1);

      //Diameter and heigth at the crown base.
      dcl.clear();
      AccumulateDown(t,dcl,AddBranchWf(),DiameterCrownBase<TS,BUD>());

      //Crown volume
      cvol = cv(t);
      //Main stem volume
      axisvol = mav(t);
      //Number of segments
      nsegment = 0;
      nsegment = AccumulateDown(t,nsegment,CountTreeSegments<TS,BUD>()); 
      qabs = 0.0;
      qabs= Accumulate(t,qabs,CollectQabs<TS,BUD>());
      treeAf = 0.0;
      treeAf = Accumulate(t,treeAf,CollectFoliageArea<TS,BUD>());
      ASeg0 = 0.0;
      ASeg0 = Accumulate(t,ASeg0,SurfaceAreaOfNewSegments());
      wood = 0.0;
      wood = Accumulate(t,wood,CollectWoodMass<TS,BUD>());

      //Collect wood mass in the main axis
      wstem = GetValue(t,LGAWstem);

      //Wood mass in branches
      wbranches = wood - wstem;

      //Sapwood mass of stem
      stem_sw = 0.0;
      Accumulate(t,stem_sw,CollectStemSapwoodMass<TS,BUD>());
      LGMdouble w_af = 0.0;
      w_af = Accumulate(t, w_af, CollectFoliageArea<TS,BUD>());

      *vdatafile[tree_n] 
	<</*1 */ left << setw(6) << setfill(' ') << iter << " "
	<</*2 */ left << setw(8) << GetValue(t,LGAH) << " " 
	<</*3 */ setw(9) << qintop1 << " " 
	<</*4 */ setw(12) << qintop1/GetValue(t,TreeQinMax) << " " 
	<</*5 */ setw(11) << GetValue(t,LGADbase) << " "
	<</*6 */ setw(11) << GetValue(t,LGADbh) << " "
	<</*7 */ setw(11) << dcl.DCrownBase() << " " 
	<</*8 */ setw(11) << dcl.HCrownBase() << " " 
	<</*9*/ setw(11)  << treePcurrent << " "      //Production before new growth
	<</*10*/ setw(11) << treeMcurrent << " "      //Respiration before new growth: foliage+sapwood+roots
	<</*11*/ setw(12) << wfaftergrowth << " " //Foliage mass after new growth
	<</*12*/ setw(14) << bs.d2l/bs.d2 << " "  //Branch means
	<</*13*/ setw(14) << bs.lsum/(double)bs.n_br  << " "
	<</*14*/ setw(11) << GetValue(t,TreeWr) << " "     //Root mass after new growth
	<</*15*/ setw(11) << GetValue(t,LGPmr)*wroot[tree_n] << " "//Root respiration before new growth
	<</*16*/ setw(11) << treeM - GetValue(t,LGPmr)*wroot[tree_n] << " "//Above ground respiration 
	<</*17*/ setw(11) << wsaftergrowth  << " "           //Sapwood mass after new growth
	<</*18*/ setw(11) << GetValue(t,LGPms)*wsapwood[tree_n] << " " //Sapwood respiration before new growth
	<</*19*/ setw(11) << GetValue(t,LGPmf)*wfoliage[tree_n] << " " //Foliage respiration before new growth
	<</*20*/ setw(11) << cvol << " "  // Crown volume
	<</*21*/ setw(11) << nsegment << " " //Number of segments
	<</*22*/ setw(11) << qabs << " "     // Qabs  
	<</*23*/ setw(11) << qabs/(GetFirmament(t).
				   diffuseBallSensor()*treeAf) << " "  //Rad eff.
	<</*24*/ setw(11) << wfoliage[tree_n] << " "  //Foliage that photosynthesized
	<</*25*/ setw(11) << treePcurrent/wfoliage[tree_n] << " "  //The ubiquitous P/Wf
	<</*26*/ setw(12) << axisvol << " "  //Main axis volume
	<</*27*/ setw(12) << ASeg0 << " "  //Main axis volume
	<</*28*/ setw(12) << wfnew << ""   //Foliage produced by elongation (iWn)
	<</*29*/ setw(12) << ws4 << " "  //Sapwood produced by elongation (iWn)
	<</*30*/ setw(12) << ws3 << " "  //Sapwood required for diameter growth (iWo)
	<</*31*/ setw(12) << wr1 << " "  //New roots required by new foliage (iWr)
	<</*32*/ setw(11) << ws4+ws3 << " " //Sapwood in growth
	<</*33*/ setw(11) << wood << " " //Wood mass in above-ground parts
	<</*34*/ setw(11) << wstem << " " //Wood mass of stem
	<</*35*/ setw(11) << wbranches << " " //Wood mass of branches
	<</*36*/ setw(11) << stem_sw << " " //Sapwood mass of stem
	<</*37*/ setw(11) << GetValue(t,LGAAsbase) << " "//Sapwood at base
	<</*38*/ setw(11) << GetValue(t,LGAAsDbh) << " " //Sapwood at D 1.3
	<</*39*/ setw(11) << dcl.ASwCrownBase() << " "       //Sapwood at crown base 
	<</*40*/ setw(11) << lambdav[0] << " " //Lambda s.t. G(L) = 0.
	<</*41*/ setw(11) << w_af << " " //Foliage area of tree
	<< endl; 
    }
  }
  
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::writeSensitivityAnalysisData(TREE& t)
  {
    if (sensitivity_analysis)
      sensitivity.printData(t);
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::writeCrownLimitData(TREE& t, int iter)
  {
    if (crown_limit_data){
      //Collect and print the two lists of foliage masses and their heights 
      //in the main axis.  Print and assess the crown base afterwards.
      CrownLimitData cld;
      AccumulateDown(t,cld,AddCrownLimitData(),
		     CollectCrownLimitData<TS,BUD>());
      //Collect diameters from the segments in the main axis
      list<TreeCompartment<TS,BUD>*> & tc_ls = GetTreeCompartmentList(GetAxis(t));
      list<double> d_ls;//list of diameters
      d_ls = accumulate(tc_ls.begin(),tc_ls.end(),d_ls,CollectSegmentDiameters());
    
      ostringstream crown_limit_file;
    
      //    crown_limit_file << "CrownLimit-" << x_coord << "-" << y_coord << "-" << iter << ".txt";
      crown_limit_file << "CrownLimit-" << GetPoint(t).getX() << "-" << GetPoint(t).getY() << "-" << iter << ".txt";

      //File is CrownLimit+iter+.txt, e.g., "CrownLimit-25.5-22.0-10.txt"
      ofstream cl_file(crown_limit_file.str().c_str());
      if (!cl_file){
	cout << "Could not open " << crown_limit_file.str() << " exit" << endl;
	exit(0);
      }
      const list<pair<double,double> >& hwf_ls = cld.WfHList();
      const list<pair<double,double> >& dh_dwf_ls = cld.dHdWfList();
      const list<pair<double,double> >& h_qabs_ls = cld.HQabsList();
      const list<pair<double,double> >& dh_dqabs_ls = cld.dHdQabsList();
      //How to do this with copy algorithm as with 'cout'??
      list<pair<double,double> >::const_iterator hwf_it = hwf_ls.begin();
      list<double>::const_iterator d_ls_it = d_ls.begin();
      list<pair<double,double> >::const_iterator dh_dwf_it = dh_dwf_ls.begin();
      list<pair<double,double> >::const_iterator h_qabs_it =  h_qabs_ls.begin();
      list<pair<double,double> >::const_iterator dh_dqabs_it = dh_dqabs_ls.begin();

      cl_file << setfill(' ') 
	      << setw(11) << "H" << " " //height
	      << setw(11) << "D" << " " //diameter
	      << setw(11) << "Wf" << " "//foliage
	      << setw(11) << "dH" << " "//growth increment 
	      << setw(11) << "dWf/dH" << " "   
	      << setw(11) << "H" << " " 
	      << setw(11) << "Qabs" << " " 
	      << setw(11) << "dH" << " " 
	      << setw(11) << "dQabs/dH" << endl; 
      while (hwf_it != hwf_ls.end()){
	cl_file << setfill(' ')
		<< setw(11) << (*hwf_it).first << " " 
		<< setw(11) << *d_ls_it << " " 
		<< setw(11) << (*hwf_it).second << " "
		<< setw(11) << (*dh_dwf_it).first << " "
		<< setw(11) << (*dh_dwf_it).second << " "
		<< setw(11) << (*h_qabs_it).first << " "
		<< setw(11) << (*h_qabs_it).second << " "
		<< setw(11) << (*dh_dqabs_it).first << " "
		<< setw(11) << (*dh_dqabs_it).second << endl;
	hwf_it++;//the heights and Wf's in the  main axis by branching points 
	dh_dwf_it++;//the same but  now dH and dWf/dH (dWf  is Wf and dH
	//the segment length)
	h_qabs_it++;//the heights and Qabs's the  in the  main axis by branching points 
	dh_dqabs_it++;//the same but now dH and dQabs/dH (dQabs is Qabs and dH
	//the segment length)
	d_ls_it++;
      }
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::writeVoxels(TREE& t)
  {
    if (interval && writevoxels){
      if (static_cast<int>(GetValue(t,LGAage)) % interval == 0){
	ostringstream voxel_space_file;
	//      voxel_space_file << "VoxelSpace-" << x_coord << "-" << y_coord << "-" << GetValue(t,LGAage) << ".txt";
	voxel_space_file << "VoxelSpace-" << GetPoint(t).getX() << "-" << GetPoint(t).getY() << "-" << GetValue(t,LGAage) << ".txt";

	vs->writeVoxelBoxesToGnuPlotFile(voxel_space_file.str());
      }
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::writeTreeToXMLFile(TREE& t, int age,int interval)const
  {
    if (verbose && xmlfile.empty()){
      cout << "writeTreeToXMLFile: No XML file given" <<endl;
    }
    if (interval && !xmlfile.empty() &&
	(age % interval == 0 )){
      Point p = GetPoint(t);
      ostringstream xml_interval;
      unsigned int  n = xmlfile.find_last_of('.');
      if ( n == 0){
	//No suffix by the user
	xml_interval << xmlfile << "-" << p.getX() <<"-"<<p.getY() << "-" << age << ".xml";
      }
      else{
	string prefix = xmlfile.substr(0,n);
	xml_interval << prefix  << "-" << p.getX() <<"-"<<p.getY() << "-" << age << ".xml";
      }
      XMLDomTreeWriter<TS,BUD> writer;
      if (verbose){
	cout << "Saving tree to "<< xml_interval.str() << " begin" <<endl; 
      }
      writer.writeTreeToXML(t,xml_interval.str());
      if (verbose){
	cout << "Saving tree end" <<endl;
      }   
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::writeFip(TREE& t,int interval)const
  {
    if (interval && !fipfile.empty()){
      if (static_cast<int>(GetValue(t,LGAage)) % interval == 0){
	ostringstream fip_interval;
	Point point = GetPoint(t);
	fip_interval << GetValue(t,LGAage) << "-" << fipfile << "-" <<GetValue(t,LGAage) <<"-"
		     << point.getX() << "-" << point.getY() << "-" << "fipfile.txt";
	if (verbose){
	  cout << "Printing vertical distribution of fip to " << fip_interval.str() << endl;
	}
	ofstream fipstream(fip_interval.str().c_str());
	int age = static_cast<int>(GetValue(t,LGAage));
	pair<vector<pair<double,int> >,double> p;
	//number of vertical divisions is tree age
	p.first.resize(age);
	//height intervals is the mean annual growth
	p.second = GetValue(t,LGAH)/GetValue(t,LGAage);
	Accumulate(t,p,CollectVerticalDistributionOfFip());
	double interval = p.second;
	double height = interval;
	//Title
	fipstream  <<left << setfill(' ')
		   << setw(11) << "Height" << setw(11) << "CumulFip" << setw(11) << "Mean fip" <<endl;
	for (unsigned int i = 0; i < p.first.size(); i++){
	  if (p.first[i].second > 0){
	    fipstream  << setw(11) << height << setw(11) <<  p.first[i].first 
		       << setw(11) << p.first[i].first/p.first[i].second <<endl;
	  }
	  else{
	    fipstream  << setw(11) << height << setw(11) <<  p.first[i].first 
		       << setw(11) << 0 <<endl;
	  }
	  height = height + interval;
	}
	if (verbose){
	  cout << "Printing fip done " <<endl;
	}
      }
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::prune()
  {
    for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
      TREE* t = vtree[k];
      LSYSTEM* l = vlsystem[k];

      l->prune(*t);
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::terminateEscapedBuds()
  {
    if (LignumForest::terminate_escaped_buds == true){
      for (unsigned int k = 0; k < (unsigned int)no_trees; k++){
	TREE* t = vtree[k];
	LSYSTEM* l = vlsystem[k];
	Point p = GetPoint(*t);
	//If a tree is not inside center stand, then it is in the border forest
	if (!terminate_buds.insideCenterStand(p)){
	  //cout << "Tree " << k << " in border stand " << p << endl; 
	  ForEach(*t,terminate_buds);
	  l->lignumToLstring(*t,1,PBDATA);
	}
      }
    }
  }
    
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::setVoxelSpaceAndBorderForest()
  {
    ///\par Voxel space bounding box
    ///Find the bounding box for the forest stand and resize voxel space.
    ///\internal
    ///\snippet{lineno} GrowthLoopI.h BoundingBox
    // [BoundingBox]
    BoundingBox bb;
    FindCfBoundingBox<TS,BUD> fb;
    for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
      bb = Accumulate(*vtree[k], bb, fb);
    }
    Point ll = bb.getMin();
    Point ur = bb.getMax();
    vs->resize(ll, ur);
    vs->reset();
    // [BoundingBox]
    ///\endinternal
    //
    ///\par Insert trees into VoxelSpace
    ///Insert trees and update voxel space  summary values.
    ///\internal
    ///\snippet{lineno} GrowthLoopI.h DumpTrees
    // [DumpTrees]
    for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
      DumpCfTree(*vs, *vtree[k], num_parts, wood_voxel);
    }
    //After dumping all trees to VoxelSpace it is necessary to evaluate
    //sum quantities (e.g. STAR_mean)
    vs->updateBoxValues();
    // [DumpTrees]
    ///\endinternal
    //
    ///\par The border forest
    ///The top is set according to voxel space (i.e. bounding box)
    ///so that the dimesions match. The dimensions of top height of the stand could
    ///also be used but it is not exactly the same as of bounding box
    ///since it considers also needles: bounding box is higher
    ///than stand top height by length of needles.
    ///\internal
    ///\snippet{lineno} GrowthLoopI.h BorderForest
    // [BorderForest]
    ///\attention Note that here only the height of BorderForest is updated. As branches
    ///\attention grow VoxelSpace extends sideways. It thus grows into BorderForest. It may be
    ///\attention necessary to correct this, depending how BorderStand.borderStandExtinction() is realised.  

    border_forest.setH(bb.getMax().getZ());
    border_forest.setHcb(stand.getMinCrownLimit());
    border_forest.setLAI(stand.getLAI());
    // [BorderForest]
    ///\endinternal
  }
 
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::calculateRadiation()
  {

    EvaluateRadiationForCfTreeSegmentInVoxelSpace
      <ScotsPineSegment,ScotsPineBud> Rad(K, vs, &border_forest, evaluate_border_forest,
					  k_border_conifer, wood_voxel, pairwise_self);

    //HUOM: true on border forest
    SetStarMean<TS,BUD> setstar(ParametricCurve(0.14));
    ResetQinQabs<TS,BUD> RQQ;
    for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
      TREE* t = vtree[k];
      ForEach(*t,setstar);
      ForEach(*t,RQQ);
  
      if(pairwise_self) {
	UnDumpScotsPineTree(*vs,*t,num_parts,wood_voxel);
      }

      ForEach(*t,Rad);

      if(pairwise_self) {
	if(k != (unsigned int)(no_trees - 1)) //Don't bother to dump the last tree,
	  //since contents of voxelspace will be erased after this
	  cout << "k " << k << endl;
	DumpCfTree(*vs, *t, num_parts, wood_voxel);
      }
    }
  }

  //end of run-voxel

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::photosynthesisRespirationTreeAging()
  {

    for (unsigned int k = 0; k < (unsigned int)no_trees; k++){
      TREE* t = vtree[k];
      photosynthesis(*t);
      respiration(*t);
      ///\par Collect data before growth
      ///
      ///Collect foliage mass and sapwood mass data before new growth.
      collectDataBeforeGrowth(*t,k);
      treeAging(*t);
      ///\par Collect sapwood mass
      ///
      ///Collect  sapwood after  senescence from  all  segments.  Collect
      ///again after  new growth excluding new  segments.  The difference
      ///of  the two  will tell how much sapwood was needed  in diameter
      ///growth.
      ws_after_senescence[k] = collectSapwoodMass(*t);
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::createNewSegments()
  {
    Pine::mode = 0;    //For  L-system
    
    for (unsigned int k = 0; k < (unsigned int)no_trees; k++){
      TREE* t = vtree[k];
      LSYSTEM* l = vlsystem[k];
      //This derive() creates the new segments, whose lengths will be iterated
      // (in allocationAndGrowth())
      l->derive();
      l->lstringToLignum(*t,1,PBDATA);
      //Pass the  qin to newly  created segments and to  the terminating
      //buds. Also set LGAip (qin/TreeQinMax) for the terminating buds
      double qin = 0.0;
      PropagateUp(*t,qin,ForwardScotsPineQin());

      //Set the needle angle for newly created segments
      ForEach(*t,SetScotsPineSegmentNeedleAngle());

      //The variable apical (now attribute of ScotsPineSegment) is set in new segments
      //before length growth by vigor index
      double o0 = 1.0;
      PropagateUp(*t,o0,SetScotsPineSegmentApical());

      TreePhysiologyVigourIndex(*t);
      //QinMax
      SetValue(*t,TreeQinMax,GetFirmament(*t).diffuseBallSensor());

      //Length of path from base of tree to each segment
      LGMdouble plength = 0.0;
      PropagateUp(*t,plength,PathLength<TS,BUD>());

      //============================================================================
      // Space colonialization
      //=============================================================================
      if (LignumForest::space0 || LignumForest::space1 || LignumForest::space2) {
	BoundingBox bbs;
	FindCfBoundingBox<TS,BUD> fbs(true); 
	//Segments with needles considered only
	bbs = Accumulate(*t, bbs, fbs);

	Point lls = bbs.getMin();
	Point urs = bbs.getMax();

	LignumForest::space_occupancy.resize(lls+Point(-0.5,-0.5,-0.5), urs+Point(0.5,0.5,0.5));
	DumpTreeOccupy<TS,BUD> dto(3);     //in 3 parts, only foliage parts
	dto.space = &LignumForest::space_occupancy;

	ForEach(*t, dto);
      }
      //==============================================================================
      // Extended Borchert-Honda calculation if it is in use
      //==============================================================================

      if(GetValue(*t, SPis_EBH) > 0.0) {
	cout << "ENTER EBH " <<endl;
	ParametricCurve lambda_fun = GetFunction(*t,SPEBHF);
	if(growthloop_is_EBH_reduction) {
	  //After age 20 EBH values for all orders change gradually from the nominal value to ebh_final_value 
	  if(L_age > 20) {
	    double v1 = lambda_fun(1.0);
	    double v2 = lambda_fun(2.0);
	    double v3 = lambda_fun(3.0);
	    double v6 = lambda_fun(6.0);
	    LGMdouble p = GetValue(*t, LGPapical);
	    v1 *= EBH_reduction_parameter;
	    v2 *= EBH_reduction_parameter;
	    v3 *= EBH_reduction_parameter;
	    v6 *= EBH_reduction_parameter;
	
	    if(EBH_reduction_parameter < 1.0) {
	      if(v1 < ebh_final_value) v1 = ebh_final_value;
	      if(v2 < ebh_final_value) v2 = ebh_final_value;
	      if(v3 < ebh_final_value) v3 = ebh_final_value;
	      if(v6 < ebh_final_value) v6 = ebh_final_value;
	    } else {
	      if(v1 > ebh_final_value) v1 = ebh_final_value;
	      if(v2 > ebh_final_value) v2 = ebh_final_value;
	      if(v3 > ebh_final_value) v3 = ebh_final_value;
	      if(v6 > ebh_final_value) v6 = ebh_final_value;
	    }

	    vector<double> lfo = lambda_fun.getVector();

	    lfo[1] = v1;
	    lfo[3] = v2;
	    lfo[5] = v3;
	    lfo[7] = v6;

	    lambda_fun = ParametricCurve(lfo);
	  }
	}
      
	EBH_basipetal_info EBHbI0, EBHbI1;
	EBHbI1 = AccumulateDown(*t, EBHbI0, EBH_basipetal(lambda_fun, LignumForest::ebh_mode) );

	EBH_acropetal_info EBHaI0(1.0, 1.0/lambda_fun(1.0), 1.0);
	PropagateUp(*t, EBHaI0, EBH_acropetal(lambda_fun) );

	MaxEBHResource_info m0, m1;
	m0.my_resource = -R_HUGE;

	m1 = AccumulateDown(*t, m0, MaxEBHResource() );

	ForEach(*t, NormalizeEBHResource(m1.my_resource) );
      }
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::radiationUseEfficiency() {
    //The maximum radiation = radiation at top of the forest =
    //ball sensor reading from the Firmament that is the same for all
    //trees is needed in calculation of radiation use efficiency of new segments.
    //It is set on the basis of shadiness experienced by their mother.
    if(growthloop_is_radiation_use_efficiency) {   //rue = 1 == no effect by default
      SetRadiationUseEfficiency<TS,BUD>  set_rue(GetFirmament(*vtree[0]).diffuseBallSensor(),
						 radiation_use_efficiency_parameter);
  
      for (unsigned int k = 0; k < (unsigned int)no_trees; k++){
	TREE* t = vtree[k];
	LGMdouble initial = 0.0;
	PropagateUp(*t,initial,set_rue);
      }
    }
  }
  
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::evaluateStandVariables() {
    stand.evaluateStandVariables(vtree, locations);
    center_stand.evaluateStandVariables(vtree, locations);
  }

  bool comp_height(pair<unsigned int,LGMdouble> p1, pair<unsigned int,LGMdouble> p2) {
    if(p1.second < p2.second)
      return true;
    else
      return false;
  }


  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::cleanUp()
  {

    //Maybe not the best place to do this but ...
    if(eero) {
      unsigned int n_trees = vtree.size();

      if(n_trees > 0) {
	vector<pair<unsigned int, LGMdouble> > sort_vector(n_trees);
	for(unsigned int i = 0; i < n_trees; i++) {
	  sort_vector[i].first = i;
	  sort_vector[i].second = GetValue(*vtree[i],LGAH);
	}

	stable_sort(sort_vector.begin(),sort_vector.end(), comp_height);

	//Take 5 percent tallest trees
	unsigned int five_per_cent = (unsigned int)(0.05 * (double)n_trees);
	if(five_per_cent < 1)
	  five_per_cent = 1;

	ofstream  f("eero_result.dat", ios_base::app);
	for(unsigned int i = n_trees - five_per_cent; i < n_trees; i++) {
	  unsigned int tree = sort_vector[i].first;
	  Point location = GetPoint(*(GetFirstTreeSegment(GetAxis(*vtree[tree]))));
	  if(center_stand.inPlot(location)) {
	    ParametricCurve SLA_fun = GetFunction(*vtree[tree], SPFSF); 
	    ParametricCurve FM_fun = GetFunction(*vtree[tree], LGMFM);
	    double treePcurrent = GetValue(*vtree[tree],TreeP);
	    double treeMcurrent = GetValue(*vtree[tree],TreeM);

	    f << location.getX() << " " << location.getY() << " " << five_per_cent 
	      << " " << GetValue(*vtree[tree],LGPpr) << " " << SLA_fun(1.0)
	      << " " << FM_fun(2.0) << " " << treePcurrent << " " << treeMcurrent << endl;

	  }
	}
	f.close();
      }
    } //  if(eero ...



    for (unsigned int i = 0; i < vlsystem.size(); i++){
      LSYSTEM* l = vlsystem[i];
      l->end();
    }

    //where this should go????
    stand_output->close();
    cstand_output->close();

  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::printSegmentQin()
  {
    for (unsigned int i = 0; i < vtree.size(); i++){
      TREE* t = vtree[i];
      //Print  Qin and  (x,y,z)  for each  segment  of age  0,1  or 2.  In
      //addition  to Qin,  this gives  you  the idea  of the  size of  the
      //branches
      if (verbose)
	cout << "Tree " << i << " Qin and point" << endl;
      ForEach(*t,PrintSegmentQin<TS,BUD>());
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::printBranchMeans()const
  {
    for (unsigned int i = 0; i < vtree.size(); i++){
      TREE* t = vtree[i];
      if (verbose)
	cout << "Tree " << i << " Branch means" << endl;
      summing bs;
      bs = Accumulate(*t, bs, Branchmeans());
      if(bs.d2 > 0.0)
	cout << "Averaged2 branch length, m: " << bs.d2l/bs.d2 << endl;
      else
	cout << "Somethin wrong in d2 branch averaging!" << endl;
    
    
      if (bs.n_br > 0)
	cout << "Average branch length, m: " << bs.lsum/(double)bs.n_br
	     << "   no. branches: " << bs.n_br << endl;
      else
	cout << "Somethin wrong in d2 branch averaging!" << endl;
    }
  }

  
  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::printTreeLocations(int iter)const
  {
    //Print the locations. Generate the file name Locations-<n>.txt where
    //n is the value of 'iter'
    ostringstream locations_file;
    locations_file << "Locations-" << iter << ".txt";
    if (verbose)
      cout << "Writing tree locations to file: " << locations_file.str() << endl;
    ofstream ofs(locations_file.str().c_str());
    for(unsigned int i = 0; i < locations.size(); i++){
      ofs << locations[i].first << " " << locations[i].second << endl;
    }
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::printVoxelObjectLocations(const string& file)const
  {
    PrintVoxelObjectLocations(*vs,file);
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::writeBranchInformation(TREE& t,const string& file)const
  {
    ForEach(t, BranchInformation(file));
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::writeProductionBalance(TREE& t,const string& file)const
  {
    ForEach(t, SegmentProductionBalance(file));
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::harvestForest(double percentage)
  {
    //Calculate number of trees to be removed
    int nremove = static_cast<int>(std::floor((percentage/100.0)*vtree.size()));
    //Collect tree positions and tree heights
    vector<pair<unsigned int,double> > vposheight;
    for (unsigned int i=0; i<vtree.size();i++){
      Tree<TS,BUD>* t = vtree[i];
      double tree_h = GetValue(*t,LGAH);
      pair<unsigned int,double> p(i,tree_h);
      vposheight.push_back(p);
    }
    //Sort the vector by tree height in ascending order
    sort(vposheight.begin(),vposheight.end(),SortByTreeHeight());
    //Collect positions to be removed
    vector<unsigned int> vremove;
    vector<pair<unsigned int,double> >::iterator vposheight_it = vposheight.begin();
    //Advance to the last position to be removed
    std::advance(vposheight_it,nremove);
    transform(vposheight.begin(),vposheight_it,std::back_inserter(vremove),CollectTreePositions());
    //Sort the positions in ascending order (default operator <)
    std::sort(vremove.begin(),vremove.end());
    //Remove trees and update data vectors
    removeHarvestedTreesAllOver(vremove);
  }

  template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::removeHarvestedTreesAllOver(const vector<unsigned int>& vremove)
  {
    cout << "Trees to be removed" << endl;
    std::copy(vremove.begin(),vremove.end(),ostream_iterator<int>(cout, " "));
    cout << endl;
    //Dead trees are removed from everywhere
    typename vector<TREE*>::iterator It = vtree.begin();
    typename vector<LSYSTEM*>::iterator Is = vlsystem.begin();
    typename vector<pair<double,double> >::iterator Il = locations.begin();
    //Vectors for data collected before new growth must be updated 
    vector<double>::iterator Iws = wsapwood.begin();
    typename vector<double>::iterator Iwf = wfoliage.begin();
    typename vector<double>::iterator Iwr = wroot.begin();
    typename vector<double>::iterator Iws_after = ws_after_senescence.begin();
    typename vector<ofstream*>::iterator If =  vdatafile.begin();
    bool also_vdatafile = false;          //vdatafile vector may be empty
    if(vdatafile.size() > 0)
      also_vdatafile = true;
    typename vector<int>::iterator In = no_h.begin();
    typename vector<double>::iterator Ih = h_prev.begin();

    typename vector<unsigned int>::const_iterator I = vremove.begin();
    unsigned int previous = 0;
    while(I != vremove.end()){
      unsigned int this_advance = *I - previous;
      advance(It, this_advance);
      advance(Is, this_advance);
      advance(Il, this_advance);
      if(also_vdatafile)
	advance(If, this_advance);
      advance(In, this_advance);
      advance(Ih, this_advance);
      //Vectors for data before new growth
      advance(Iws,this_advance);
      advance(Iwf,this_advance);
      advance(Iwr,this_advance);
      advance(Iws_after,this_advance);
	
      delete *It;   //locations were not created by new
      delete *Is;
      if(also_vdatafile)
	delete *If;


      cout << "Dead tree at location " << Il->first << " " << Il->second << " was deleted in year "
	   << year << endl;

      vtree.erase(It);
      vlsystem.erase(Is);
      locations.erase(Il);
      //Vectors  for data  before new  growth must  also be  udated to
      //maintain the integrity of  the positions (indices) of trees to
      //these vectors
      wsapwood.erase(Iws);
      wfoliage.erase(Iwf);
      wroot.erase(Iwr);
      ws_after_senescence.erase(Iws_after);
      if(also_vdatafile)
	vdatafile.erase(If);
      no_h.erase(In);
      h_prev.erase(Ih);

      previous = *I + 1;

      I++;
      no_trees--;
    }
  }
}//End namespace LignumForest
#endif
