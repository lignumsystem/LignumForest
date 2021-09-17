#ifndef GROWTHLOOPI_H
#define GROWTHLOOPI_H


// Tama GrowthLoopI tiedosto ajaa kokonaan voxelspaessa, muutokset
// merkattu: run-voxel

namespace Pine{

  extern int mode;  //For  L-system; change of mode disabled there (/LignumForest/pine-em98.L)

}
using namespace Pine;


extern int ran3_seed;
extern double H_0_ini, H_var_ini;         //For variation of initial heights and
extern int n_buds_ini_min, n_buds_ini_max;  // and number of buds (in .L file)

extern double L_age;

extern double rel_bud;                   //Variation in no. buds
extern bool bud_variation;             //If bud variation is on
extern double branch_angle;


/* template <class TS, class BUD> */
/*   class EvaluateRadiationForCfTreeSegment_3; */




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
  
  stand_output->close();
  delete stand_output;
  
}

template<class TREE, class TS, class BUD, class LSYSTEM>
void GrowthLoop<TREE,TS,BUD,LSYSTEM>::usage()const
{
  cout << "Usage:  ./lig-forest -iter <value>  -metafile <file>  -voxelspace <file>" <<endl;
  cout << "[-numParts <parts>]  [-treeDist <dist>] [-hw <hw_start>] [-viz]" <<endl;
  cout << "[-toFile <filename>] [-xml <filename>] [-writeVoxels] [-sensitivity <filename>] " <<endl;
  cout << "[-fipdistrib <filename>] [-writeInterval interval]" << endl;
  cout << "[-seed <num>] [-increaseXi <value>] [-targetTree <num>] " <<endl;
  cout << "[-treeFile <filename>] [-generateLocations  <num>] [-woodVoxel] [-treeLocations <file>]" << endl;
  cout << "[-writeOutput] [-verbose] [-bracketVerbose] [-noBorderForest] [-seed <value>] [-kBorderConifer <value>]"  << endl;
  cout << "[-H_0_ini <value>] [-H_var_ini <value>] [-n_buds_ini_min <num>] [-n_buds_ini_max <vlaue>]" << endl;
  cout << "[-p0Var <value>] [-segLenVar <value>] [-pairwiseSelf] [-budVariation <value>] [-eero]" << endl;
  cout << "[-gFunVar <value>] [-branchAngleVar <value>]" << endl;

  cout << endl;
  cout << "-generateLocations <num>  In this case <num> trees will be generated to random locations. If this" << endl;
  cout << "          is not on, tree locations will be read from file Treelocations.txt. This file can be changed" << endl;
  cout << "          by -treeLocations <file>. If location file is not found program stops." << endl;
  cout << "-woodVoxel                If woody parts are dumped to voxels (default = true)" << endl;
  cout << "-treeDist <dist>          Minimum distance between two trees (default = 0), works only with -generateLocations." << endl;  
  cout << "-numParts <parts>         Segments can be dumped to voxels by parts (i.e. they may belong to different voxels," << endl;
  cout << "                          default = 1)" << endl;
  cout << "-hw <hw_start>            Starting year of formation of sapwood. Default = 15 years" << endl;
  cout << "-increaseXi               If parameter ksi (fraction of heartwood in new segments) increases, value = starting year (default = 15)"
       << endl;
  cout << "                          as increase = 0.004/year after year 15 up to value 0.85 (as in FPB 35: 964-975 2008 article)"
       << endl;
  cout << "-targetTree <num>         Any one of the trees can be identified as target tree (default = 0)" << endl;
  cout << "-writeOutput                Most of the things are written to their respctive file at -writeInterval interval (default false)" << endl;
    cout << "-verbose                  Output of progress of run if set (default = not set)." << endl;
  cout << "-bracketVerbose             If set, iteration information is printed out in allocation (default = not set)." << endl;
  cout << "-dumpSelf                   If the subject tree is dumped to vox-space (default no)" << endl;
  cout << "-noBorderForest             No border forest around the stand (default = there is border forest)"
       << endl;
  cout << "-seed <value>               seed for random number generator." << endl;
  cout << "-kBorderConifer <value>       Extinction coefficient for conifer foliage in border forest (default = 0.14)" << endl;
  cout << "-H_0_ini, -H_var_ini        For variation of initial heights (defaults = 0.3 and 0.0)" << endl;
  cout << "-n_buds_ini_min, -n_buds_ini_max  For variation of initial number of buds (defaults = 4 and 4)" << endl;
  cout << "-p0Var <value>                Random variation of p0 +- max <value> per cent from the value in Tree.txt" << endl;
  cout << "-segLenVar <value>         Random variation of length of new segments around Lnew, per cent" << endl;
  cout << "-pairwiseSelf      Pairwise radiation calculation for the tree itself." << endl;

  cout << endl;

}



template<class TREE, class TS, class BUD, class LSYSTEM>
void GrowthLoop<TREE,TS,BUD,LSYSTEM>::checkCommandLine(int argc, char** argv)const
{
  //At least three  mandatory arguments required 
  if (argc < 4){
    cout << "Three mandatory arguments are required!" << endl << endl;
    usage();
    exit(0);
  }
  else if (CheckCommandLine(argc,argv,"-iter") == false){
    cout << "Mandatory -iter <num> option missing" << endl;
    exit(0);
  }
  else if (CheckCommandLine(argc,argv,"-metafile") == false){
    cout << "Mandatory -metafile <MetaFile.txt> option missing" << endl;
    exit(0);
  }
  else if (CheckCommandLine(argc,argv,"-voxelspace") == false){
    cout << "Mandatory -voxelspace <VoxelSpace.txt> option missing" << endl;
    exit(0);
  }
  else if (verbose){
    cout << "Command line O.K." <<endl;
  } 
}

template<class TREE, class TS, class BUD, class LSYSTEM>
void GrowthLoop<TREE,TS,BUD,LSYSTEM>::parseCommandLine(int argc, char** argv)
{
  //first verbose
  verbose = false;
  if (CheckCommandLine(argc,argv,"-verbose")){
    verbose = true;
  }

  if (verbose){
    cout << "parseCommandLine begin" <<endl;
  }

  checkCommandLine(argc,argv);

  //Mandatory arguments
  string clarg;
  if (ParseCommandLine(argc,argv,"-iter", clarg)){
    iterations = atoi(clarg.c_str());
  }

  clarg.clear();
  if (ParseCommandLine(argc,argv,"-metafile", clarg)){
    metafile = clarg;
  }

  clarg.clear();
  if (ParseCommandLine(argc,argv,"-voxelspace", clarg)){
    voxelfile = clarg;
  }

  //End of mandatory arguments 

//  clarg.clear();
//  if (ParseCommandLine(argc,argv,"-treeDist", clarg))
//    tree_distance = atof(clarg.c_str());


  //Output to data file. This is  the base name each tree will add its
  //coordinates to make unique files
  clarg.clear();
  to_file = ParseCommandLine(argc,argv,"-toFile", clarg);
  if (to_file){
    datafile = clarg;
  }

  cstand_file = "cstand-values.dat"; 
  clarg.clear();
  ParseCommandLine(argc,argv,"-cstandFile", clarg);
    cstand_file = clarg;


stand_file = "stand-values.dat"; 
clarg.clear();
  if(ParseCommandLine(argc,argv,"-standFile", clarg))
    stand_file = clarg;
 //sensitivity analysis file
  clarg.clear();
  sensitivity_analysis = ParseCommandLine(argc,argv,"-sensitivity", clarg);
  if (sensitivity_analysis){
    sensitivity.printHeader(clarg);
  }

  //Write  crown limit data
  crown_limit_data =  CheckCommandLine(argc,argv,"-crowmLimitData");

  //XML file where the tree can be saved and restored from
  clarg.clear();
  ParseCommandLine(argc,argv,"-xml", clarg);
  xmlfile = clarg;

  //The vertical distribution of fip from segments
  clarg.clear();
  ParseCommandLine(argc,argv,"-fipdistrib", clarg);
  fipfile = clarg;

  //Write voxels, the file will be named as 'VoxelSpace-x-y-age.txt'
  if (CheckCommandLine(argc,argv,"-writeVoxels")){
    writevoxels = true;
  }

  //Number of segment parts used to assess Qabs
  num_parts = 1;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-numParts", clarg))
    num_parts = atoi(clarg.c_str());
  
  wood_voxel = true;
  if (CheckCommandLine(argc,argv,"-noWoodVoxel"))
    wood_voxel = false;

  //Write interval in output
  interval = 1;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-writeInterval", clarg))
    interval = atoi(clarg.c_str());

  write_output = false;
  if (CheckCommandLine(argc,argv,"-writeOutput"))
    write_output = true;

		     
  //Increase LGPxi value

  clarg.clear();
  if (ParseCommandLine(argc,argv,"-increaseXi", clarg)){
    increase_xi = true;
    xi_start = atoi(clarg.c_str());
  }

  //The start of heartwood build up
  clarg.clear();
  hw_start = 15;
  if (ParseCommandLine(argc,argv,"-hw", clarg))
    hw_start = atoi(clarg.c_str());


  //Initialize ran3
  //clarg.clear();
  //if (ParseCommandLine(argc,argv,"-seed", clarg)){
  //  int s = atoi(clarg.c_str());
  //  initRan3(s);
//  }

  generate_locations = false;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-generateLocations", clarg)) {
    generate_locations = true;
    no_trees = atoi(clarg.c_str());
  }
  else {
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-treeLocations", clarg)) 
      location_file = clarg;
    else
      location_file = "Treelocations.txt";
  }

  tree_distance = 0.0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-treeDist", clarg))
    tree_distance = atof(clarg.c_str());

  target_tree = 0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-targetTree", clarg))
    target_tree = (unsigned int)atoi(clarg.c_str());


  dump_self = false;
  if (CheckCommandLine(argc,argv,"-dumpSelf"))
    dump_self = true;


  bracket_verbose = false;
  if (CheckCommandLine(argc,argv,"-bracketVerbose"))
    bracket_verbose = true;

  evaluate_border_forest = true;
  if (CheckCommandLine(argc,argv,"-noBorderForest"))
    evaluate_border_forest = false;

  ran3_seed = -123231;
  int s_ini;
  if (ParseCommandLine(argc,argv,"-seed", clarg)){
    if (clarg.length() > 0){
      s_ini = atoi(clarg.c_str());
      ran3_seed = -abs(s_ini);
    }
  }

  k_border_conifer = 0.14;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-kBorderConifer", clarg))
    k_border_conifer = atof(clarg.c_str());

  H_0_ini = 0.3;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-H_0_ini", clarg))
    H_0_ini = atof(clarg.c_str());

  H_var_ini = 0.0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-H_var_ini", clarg))
    H_var_ini = atof(clarg.c_str());

  n_buds_ini_min = 4;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-n_buds_ini_min", clarg))
    n_buds_ini_min  = atoi(clarg.c_str());

  n_buds_ini_max = 4;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-n_buds_ini_max", clarg))
    n_buds_ini_max = atoi(clarg.c_str());

  p0_var = 0.0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-p0Var", clarg))
    p0_var = atof(clarg.c_str());

  seg_len_var = 0.0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-segLenVar", clarg))
    seg_len_var = atof(clarg.c_str());
  seg_len_var /= 100.0;            //per cent -> [0,1]

  pairwise_self = false;
  if (CheckCommandLine(argc,argv,"-pairwiseSelf")) {
     pairwise_self = true;
  }

  rel_bud = 0.0;
  bud_variation = false;
 clarg.clear();
 if (ParseCommandLine(argc,argv,"-budVariation", clarg)) {
   bud_variation = true;
   rel_bud = atof(clarg.c_str())/100.0;
 }

 eero = false;
  if (CheckCommandLine(argc,argv,"-eero")) {
     eero = true;
  }

  g_fun_var = 0.0;
  g_fun_varies = false;
 clarg.clear();
 if (ParseCommandLine(argc,argv,"-gFunVar", clarg)) {
   g_fun_var = atof(clarg.c_str())*PI_VALUE/100.0;
   g_fun_varies = true;
 }

 random_branch_angle = false;
 if (ParseCommandLine(argc,argv,"-branchAngleVar", clarg)) {
   random_branch_angle = true;
   ba_variation = atof(clarg.c_str())*PI_VALUE/100.0;
 }


  if (verbose){
    cout << "parseCommandLine end" <<endl;
    printVariables();
  }   
}




template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::resolveCommandLineAttributes()
{
  if(pairwise_self) {
    dump_self = false;
  }

  if(eero) {
    p0_var = 0.0;
    bud_variation = false;
    seg_len_var = 0.0;
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

//==============================================================================================
//Create the  trees.
//===============================================================================================
template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::createTrees()
{
  for (int i = 0; i < no_trees; i++){
    pair<double,double> p = locations[i];
    LSYSTEM *l = new LSYSTEM();
    TREE* t = new TREE(Point(p.first,p.second,0.0),PositionVector(0,0,1),
		       "sf.fun","fapical.fun","fgo.fun",
		       "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
		       "flr.fun");

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
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::initializeTrees()
{
  LGMVERBOSE vrb=QUIET;
  if (verbose)
    vrb = VERBOSE;
  InitializeTree<TS,BUD> init(metafile,vrb);
  for (unsigned int i=0; i < vtree.size(); i++){
    SetValue(*vtree[i],SPHwStart,(double)hw_start);
    init.initialize(*vtree[i]);
    SetValue(*vtree[i],TreeQinMax,GetFirmament(*vtree[i]).diffuseBallSensor());
    SetValue(*vtree[i], LGPlen_random, seg_len_var);
  }

  //Random variation in the photosynthetic parameter

  for (unsigned int i=0; i < vtree.size(); i++){
    LGMdouble p0 = GetValue(*vtree[i],LGPpr);
    p0 *= 1.0 + ((ran3(&ran3_seed) - 0.5)/0.5)*p0_var/100.0;
    SetValue(*vtree[i],LGPpr,p0);
  }
 
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

      //Light-use efficiency
      LGMdouble LUE = GetValue(*vtree[i],LGPpr);
      LGMdouble LUE_change = 1.0 + ((ran3(&ran3_seed)-0.5)/0.5) * LUE_factor;
      if(LUE_change < 0.0)
	LUE_change = 0.0;

      SetValue(*vtree[i],LGPpr,LUE_change*LUE);
      LGMdouble mf = GetValue(*vtree[i],LGPmf);
      SetValue(*vtree[i],LGPmf,LUE_change*mf);      //LUE and respiration correlated

      //SLA - assumes that the function is defined at argument values 0.0, 1.0, and 1.1
      LGMdouble SLA_change = 1.0 + ((ran3(&ran3_seed)-0.5)/0.5) * SLA_factor;
      while(SLA_change < 0.0) {
	 SLA_change = 1.0 + ((ran3(&ran3_seed)-0.5)/0.5) * SLA_factor;
      }
/*       if(SLA_change < 0.0) */
/* 	SLA_change = 0.0; */
      ParametricCurve SLA_fun = GetFunction(*vtree[i], SPFSF);

      ostringstream to_new_SLA_fun;
      to_new_SLA_fun << 0.0 << " " << SLA_change*SLA_fun(0.0) << " " << 1.0 << " "
		     << SLA_change*SLA_fun(1.0) << " " << 1.1 << " " << SLA_change*SLA_fun(1.1);

      int dummy = 0;
      ParametricCurve new_SLA_fun(to_new_SLA_fun.str(), dummy);
      SetFunction(*vtree[i], new_SLA_fun, SPFSF);

      //Foliage mortality (=foliage longevity) - assumes that the function has been
      //defined at values 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, and 6.0
      //Foliage longevity changes are obtained by stretching the argument values 

      LGMdouble FMc = 1.0 + ran3(&ran3_seed)*fol_age_factor;
      //      int FMc = (int)FM_change;
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
    }
  }   //eero

  if(g_fun_varies) {
    for (unsigned int i=0; i < vtree.size(); i++){

      ParametricCurve G_fun = GetFunction(*vtree[i],SPFGO);
      //Assumes that function has been defined for orders 0, ..., 7, orders 0,1 are
      //not changed
      LGMdouble G_c = 1.0 + ((ran3(&ran3_seed)-0.5)/0.5) * g_fun_var;
      double v2 = G_c*G_fun(2.0);
      double v3 = G_c*G_fun(3.0);
      double v4 = G_c*G_fun(4.0);
      double v5 = G_c*G_fun(5.0);
      double v6 = G_c*G_fun(6.0);
      double v7 = G_c*G_fun(7.0);
      if(v2 > 1.0) v2 = 1.0;
      if(v3 > 1.0) v3 = 1.0;
      if(v4 > 1.0) v4 = 1.0;
      if(v5 > 1.0) v5 = 1.0;
      if(v6 > 1.0) v6 = 1.0;
      if(v7 > 1.0) v7 = 1.0;
      ostringstream to_new_G_fun;
      to_new_G_fun << 0.0 << " " << G_fun(0.0) << " "
		   << 1.0 << " " << G_fun(1.0) << " "
		   << 2.0 << " " << v2 << " "
		   << 3.0 << " " << v3 << " "
		   << 4.0 << " " << v4 << " "
		   << 5.0 << " " << v5 << " "
		   << 6.0 << " " << v6 << " "
		   << 7.0 << " " << v7;
      int dummy = 0;
      ParametricCurve new_G_fun(to_new_G_fun.str(), dummy);
      SetFunction(*vtree[i], new_G_fun, SPFGO);
    }
  }     //if(g_fun_varies)


  //must set branch_angle to trees any way
  branch_angle = 45.0 * PI_VALUE / 180.0;

  for (unsigned int i=0; i < vtree.size(); i++){
    if(random_branch_angle) {
      branch_angle *= 1.0 + ((ran3(&ran3_seed)-0.5)/0.5) * ba_variation; 
    }
    vtree[i]->setBranchAngle(branch_angle);
  }  

}

template<class TREE, class TS,class BUD, class LSYSTEM>
void GrowthLoop<TREE, TS,BUD,LSYSTEM>::initializeVoxelSpace()
{
  ifstream vf(voxelfile.c_str());
  LGMdouble vx,vy,vz; vx = vy = vz = 0.0;
  LGMdouble s1,s2,s3,b1,b2; s1 = s2 = s3 = 0.0; b1=b2=0.0;
  vf >> vx >> vy >> vz >> s1 >> s2 >> s3 >> b1 >> b2;
  //vx,vy,vz define the voxel space  dimensions, s1, s2, s3 define the
  //voxel box dimensions
  if (verbose){
    cout << "Voxel Space: " << vx << " " << vy << " " << vz << " " 
	 << s1 << " " << s2 << " " << s3 << " " << b1 << " " << b2 <<endl;
  }
  vs = new VoxelSpace(Point(0,0,0),Point(vx,vy,vz),
		      s1,s2,s3,
		      static_cast<int>(vx/s1),static_cast<int>(vy/s2),static_cast<int>(vz/s3),
		      GetFirmament(*vtree[0]));
}

template<class TREE, class TS,class BUD, class LSYSTEM>
void GrowthLoop<TREE, TS,BUD,LSYSTEM>::initializeFunctions()
{  
  K.install("K.fun");
  stems_ha.install("stemsha.fun");
  fdensity_ha.install("fdensity.fun");
  if (verbose){
    cout << "K() O.K: " << K.ok() << " stems_ha() O.K: " << stems_ha.ok() 
	 << " fdensity() O.K: " << fdensity_ha.ok() << endl;
  }
}
//================================================================================
//Generate tree locations, or read them from a file
//Establish also stand corners with this information
//They are set in StandDescriptor and BorderForest
//================================================================================
template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::setTreeLocations()
{  
  if(generate_locations) {
    //In this case the take the plot dimensions from voxelspace file: Read the voxel space file
    ifstream vf(voxelfile.c_str());
    LGMdouble vx,vy,vz; vx = vy = vz = 0;
    LGMdouble s1,s2,s3,b1,b2; s1 = s2 = s3 = 0.0; b1=b2=0.0;
    vf >> vx >> vy >> vz >> s1 >> s2 >> s3 >> b1 >> b2;
    
    //Corners of the stand
    Point l(0.0, 0.0, 0.0);
    Point r(vx, vy, 0.0);
    stand.setLlCorner(l);
    stand.setUrCorner(r);
    stand.evaluateArea();

    //Set the same for center_stand that is used to evaluate stand values
    //without border effect
    Point cl(b1,b2,0.0);
    Point cr(vx-b1, vy-b2, 0.0);
    center_stand.setLlCorner(cl);
    center_stand.setUrCorner(cr);
    center_stand.evaluateArea();

    //Border Forest
    border_forest.setCornerL(l);
    border_forest.setCornerR(r);

    //ForestGap is here only for consistency with use of GenerateLocations in Lig-Crobas
    //    gap_radius = 0.0;
    //        ForestGap gap(pair<double,double>(x_coord,y_coord),gap_radius);
        ForestGap gap(pair<double,double>(0.0,0.0),0.0);

    //number of trees may decrease due to hard core
    int no_trees_0 = no_trees;
        GenerateLocations(no_trees,0.0,0.0,vx,vy,tree_distance,gap,locations);
    //Insert the target tree location
	//    locations.insert(locations.begin(), pair<double,double>(x_coord,y_coord));
    if (verbose){
      cout << "Number of trees" << locations.size() <<endl 
	   << " Density/ha wanted: " << (double)no_trees_0/(vx*vy/10000.0)
           << " Density/ha created: " << (double)no_trees/(vx*vy/10000.0) <<endl;
      cout << " Minimum tree distance: " << tree_distance <<endl; 
    }
  } //  if(generate_  ...)
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
      double xii = GetValue(*t, LGPxi);
      xii += 0.1/25.0;   
      if(xii > 0.85) xii=0.85;
      SetValue(*t,LGPxi, xii);
      if (verbose && first) {
	cout << "Increasing LGPxi from " << GetValue(*t,LGPxi) << " to " << xii <<endl;
	first = false;
      }
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
  //To print  consistent P and  respirations we must  collect masses
  //now before senescense and new growth
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

template<class TREE, class TS,class BUD, class LSYSTEM>
void GrowthLoop<TREE, TS,BUD,LSYSTEM>::treeAging(TREE& t)
{
  //TreeAging takes care of senescence in segments (above ground part)
    ForEach(t,TreeAging<TS,BUD>()); 
    //Root mortality
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

//========================================================================
// allocation() does the allocation, that is, finding value of lambda that
// makes demand to match available resources: P _ M = use in growth
// with the aid of iteration.
// In the case P - M < 0 or iteration cannot find solution, the tree
// can be considered to be dead, and  allocation() returns false,
// otherwise it returns true.
//========================================================================

template<class TREE, class TS,class BUD, class LSYSTEM>
bool GrowthLoop<TREE, TS,BUD,LSYSTEM>::allocation(TREE& t, bool verbose)

{
  //Allocate    net    photosynthesis    
  //Testing the implementation where the sapwood area is passed down
  //as such  between segments that are  in the same  axis.  Only the
  //segments  of higher gravelius  order require  less sapwood  in a
  //branching point.

  DiameterGrowthData data;
  LGMGrowthAllocator2<TS,BUD,SetScotsPineSegmentLength,
      PartialSapwoodAreaDown,ScotsPineDiameterGrowth2,DiameterGrowthData>
    G(t,data,PartialSapwoodAreaDown(GetFunction(t,SPSD))); 
 
    treeP = G.getP();
    treeM = G.getM();

   if(treeP - treeM < 0.0)
    return false;

  try{
    Bisection(0.0,10.0,G,0.01,verbose); //10 grams (C) accuracy 
    lambda = G.getL();
  }
  //G will throw an exception if P < M
  catch(TreeGrowthAllocatorException e){
    cout << "P < M " << e.getP() << " " << e.getM() <<endl;
    return false;
  }
  catch(BisectionBracketException e){
    cout << "Could not bracket " << e.getFa() << " " << e.getFb() << " "  << e.getFbl()  <<endl;
    cout << e.getA() << " "  << e.getB() << " " << e.getBl() <<endl;
    return false;
  }

  return true;
}


//===================================================================
// Allocation of photosynthetic production to growth of new segments and
// expansion of existing ones. The tree structure is also updated and
// new buds are created.
// If iterative allocation does not succeed in allocation() (probably
// P - M < 0.0) the tree is considered dead and removed from the tree
// list (and no_trees is updated)
//===================================================================
template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::allocationAndGrowth()
{
  mode = 1;     //For  L-system; change of mode disabled there (/LignumForest/pine-em98.L)

  list<unsigned int> dead_trees;
  dead_trees.clear();
 
  for (unsigned int k = 0; k < (unsigned int)no_trees; k++){
    TREE* t = vtree[k];
    LSYSTEM* l = vlsystem[k];

    if(!allocation(*t,bracket_verbose))
 
      dead_trees.push_back(k);       //iteration failed, this tree is dead
    else
      if(GetValue(*t,LGAH) - h_prev[(int)k] < 0.001) {
	cout << L_age << " " << k << " " << no_h[(int)k] << endl;
	no_h[(int)k] += 1;
	if(no_h[(int)k] >= 3)
	   dead_trees.push_back(k);
      }

    //Calculate  the  LGAsf  for   newly  created  segments,  sf  in  P
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

    //Create new buds by making derive again (done in createNewSegments()). This
    // works since l.derive() makes new buds with every second call
    branch_angle = t->getBranchAngle();        //this global variable goes to pine-em98.L
    l->derive();
    l->lstringToLignum(*t,1,PBDATA);
  }

  //Now, dead trees are removed from everywhere
  typename vector<TREE*>::iterator It = vtree.begin();
  typename vector<LSYSTEM*>::iterator Is = vlsystem.begin();
  vector<pair<double,double> >::iterator Il = locations.begin();
  //Vectors for data collected before new growth must be updated 
  vector<double>::iterator Iws = wsapwood.begin();
  vector<double>::iterator Iwf = wfoliage.begin();
  vector<double>::iterator Iwr = wroot.begin();
  vector<double>::iterator Iws_after = ws_after_senescence.begin();
  vector<ofstream*>::iterator If =  vdatafile.begin();
  bool also_vdatafile = false;          //vdatafile vector may be empty
  if(vdatafile.size() > 0)
    also_vdatafile = true;
  typename vector<int>::iterator In = no_h.begin();
  typename vector<double>::iterator Ih = h_prev.begin();

  list<unsigned int>::iterator I = dead_trees.begin();
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

//This method will print various data of a tree into a given file. If you call this method
//after GrowthLoop::allocationAndGrowth() then the tree has newly created segments. The GrowthLoop::treeP
//and GrowthLoop::treeM 
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
      <</*40*/ setw(11) << lambda << " " //Lambda s.t. G(L) = 0.
      <</*41*/ setw(11) << w_af << " " //Foliage area of tree
      << endl; 



    if(pairwise_self) {
      if(iter == 9){
	Point p = GetPoint(*GetFirstTreeSegment(GetAxis(t)));
	LGMdouble xx = p.getX();
	LGMdouble yy = p.getY();
	if((abs(xx-12.2775) < 0.0001) && (abs(yy-7.57141)<0.0001))
	  ForEach(t, SegmentProductionBalance("phprod-10.dat"));
      }

      if(iter == 14){
	Point p = GetPoint(*GetFirstTreeSegment(GetAxis(t)));
	LGMdouble xx = p.getX();
	LGMdouble yy = p.getY();
	if((abs(xx-12.2775) < 0.0001) && (abs(yy-7.57141)<0.0001))
	  ForEach(t, SegmentProductionBalance("phprod-15.dat"));
      }
      if(iter == 19){
	Point p = GetPoint(*GetFirstTreeSegment(GetAxis(t)));
	LGMdouble xx = p.getX();
	LGMdouble yy = p.getY();
	if((abs(xx-12.2775) < 0.0001) && (abs(yy-7.57141)<0.0001))
	  ForEach(t, SegmentProductionBalance("phprod-20.dat"));
      }


    }



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
      cout << "Could not open " << crown_limit_file << " exit" << endl;
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

//Write a tree into xml file.  The file name comprises of the basename
//(=xmlfile), the position of the tree (x, y), and its age.
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

//Vertical distrubution of fip
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
//===================================================================
// Prepare for radiation calculations:
// Set BoundingBox and VoxelSpace for trees in the stand.
// Dump also the foliage (except the first tree) into voxels and
// set BorderForest
//===================================================================
template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::setVoxelSpaceAndBorderForest()
{
  BoundingBox bb;
  FindCfBoundingBox<TS,BUD> fb;
  for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
    bb = Accumulate(*vtree[k], bb, fb);
  }

  Point ll = bb.getMin();
  Point ur = bb.getMax();

  vs->resize(ll, ur);

  vs->reset();

  //Note: not first tree, since it will be the first to be calculated
  //and its foliage won't be in the voxelspace
 // unless dump_self is set

  //  for (unsigned int k = 1; k < (unsigned int)no_trees; k++)
  if(dump_self)
      DumpCfTree(*vs, *vtree[0], num_parts, wood_voxel);
  for (unsigned int k = 1; k < (unsigned int)no_trees; k++)
    DumpCfTree(*vs, *vtree[k], num_parts, wood_voxel);

  if((int)GetValue(*vtree[0],LGAage) == 10){
    LGMdouble Hmin, Hmax;
    int nz;
    vector<pair<double,double> > NAD;
    vs->evaluateVerticalNeedleAreaDensity(Hmax, Hmin, nz, NAD);
    cout << "DISTN DISTN" << endl;
    cout << "nz Hmin Hmax " << nz << " " << Hmin << " " << Hmax << endl;
    for(int i = 0; i < nz; i++) {
      cout << NAD[i].first << " " << NAD[i].second << endl;
    }
  }

  if((int)GetValue(*vtree[0],LGAage) == 15){
    LGMdouble Hmin, Hmax;
    int nz;
    vector<pair<double,double> > NAD;
    vs->evaluateVerticalNeedleAreaDensity(Hmax, Hmin, nz, NAD);
    cout << "DISTN DISTN" << endl;
    cout << "nz Hmin Hmax " << nz << " " << Hmin << " " << Hmax << endl;
    for(int i = 0; i < nz; i++) {
      cout << NAD[i].first << " " << NAD[i].second << endl;
    }
  }
  if((int)GetValue(*vtree[0],LGAage) == 20){
    LGMdouble Hmin, Hmax;
    int nz;
    vector<pair<double,double> > NAD;
    vs->evaluateVerticalNeedleAreaDensity(Hmax, Hmin, nz, NAD);
    cout << "DISTN DISTN" << endl;
    cout << "nz Hmin Hmax " << nz << " " << Hmin << " " << Hmax << endl;
    for(int i = 0; i < nz; i++) {
      cout << NAD[i].first << " " << NAD[i].second << endl;
    }
  }


  //After dumping all trees to VoxelSpace it is necessary to evaluate
  //sum quantities (e.g. STAR_mean)
  vs->updateBoxValues();


  //The border forest top is set according to voxelspace
  // (i.e. bounding box)
  // so that the dimesions match. The dimensions of top height of the stand could
  // also be used but it is not exactly the same as of bounding box
  // since it considers also needles: bounding box is higher
  // than stand top height by length of needles. 
  border_forest.setH(bb.getMax().getZ());
  border_forest.setHcb(stand.getMinCrownLimit());
  border_forest.setLAI(stand.getLAI());
}

//===================================================================
// calculateRadiation UnDumps a tree (except the first tree, since
// it was not dumped in setVoxelspace()) and calculates radiation
// 1) by pairwise for itself and 2) through voxelspace and
// 3) borderforest outside itself
// and then dumps it back to voxelspace
//===================================================================
template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::calculateRadiation()
{

  LGMdouble a = 1.0, b = 1.0;
  EvaluateRadiationForCfTreeSegment_3
    <ScotsPineSegment,ScotsPineBud> Rad(K, vs, &border_forest, evaluate_border_forest, a, b, dump_self,
					k_border_conifer, false, false, false, false, false);
/*    EvaluateRadiationForCfTreeSegment_3<ScotsPineSegment,ScotsPineBud> */
/*      Rad3(K, vs, &border_forest, false, a, b, virittely_dump, k_border_conifer, */
/*          box_dir_effect, wood_voxel, correct_star, constant_star,calculateDirectionalStar); */


  //HUOM: true on border forest

  SetStarMean<TS,BUD> setstar(ParametricCurve(0.14));
  ResetQinQabs<TS,BUD> RQQ;

  for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
    TREE* t = vtree[k];
    ForEach(*t,setstar);
    ForEach(*t,RQQ);
    
    if(pairwise_self) {
      if(k != 0)       //First tree was not dumped in
	// setting up voxelspace in setVoxelSpaceAndBorderForest
	UnDumpScotsPineTree(*vs,*vtree[k],num_parts,wood_voxel);
    }

    ForEach(*t,Rad);

    if(pairwise_self) {
      if(k != (unsigned int)(no_trees - 1)) //Don't bother to dump the last tree,
	//since contents of voxelspace will be erased after this
	DumpCfTree(*vs, *t, num_parts, wood_voxel);
    }

  }

}

//end of run-voxel

//===================================================================
// Photosynthesis and respiration & collection of some mass variables
//===================================================================
template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::photosynthesisAndRespiration()
{

  for (unsigned int k = 0; k < (unsigned int)no_trees; k++){
    TREE* t = vtree[k];
    photosynthesis(*t);
    respiration(*t);
    collectDataBeforeGrowth(*t,k);
    treeAging(*t);
    //Collect  sapwood after  senescence from  all  segments.  Collect
    //again after  new growth excluding new  segments.  The difference
    //of  the two  will tell  how much  sapwood was  need  in diameter
    //growth
    ws_after_senescence[k] = collectSapwoodMass(*t);

  }
}

//===================================================================
// Create new segments and set some some variables
// NOTE that the sizes of these new segments will be iterated later
//===================================================================
template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::createNewSegments()
{
  mode = 0;    //For  L-system; change of mode disabled there (/LignumForest/pine-em98.L)
     
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

    TreePhysiologyVigourIndex(*t);

    //QinMax
    SetValue(*t,TreeQinMax,GetFirmament(*t).diffuseBallSensor());

    LGMdouble plength = 0.0;
    PropagateUp(*t,plength,PathLength<TS,BUD>());

    double alku = 1.0;    //= Gravelius order of main axis
    PropagateUp(*t,alku,SetSapwoodDemandAtJunction());
  }

}

//=====================================================================
// Evaluate stand variables of stand and center_stand
//=====================================================================

template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::evaluateStandVariables() {

  stand.evaluateStandVariables(vtree, locations);
  center_stand.evaluateStandVariables(vtree, locations);
}




//===================================================================
// Output all
//===================================================================
template<class TREE, class TS,class BUD, class LSYSTEM>
  void GrowthLoop<TREE, TS,BUD,LSYSTEM>::output()
{
  if(write_output) {
      *stand_output << year << " ";
      stand.writeOutput(*stand_output);

      *cstand_output << year << " ";
      center_stand.writeOutput(*cstand_output);


    for (unsigned int k = 0; k < (unsigned int)no_trees; k++){
      TREE* t = vtree[k];
      //Annual (time  step) output: tree, its position  in the vector
      //and iteration
      writeOutput(*t,k,year);
      writeCrownLimitData(*t,year);
      //      writeVoxels(*t);
      writeTreeToXMLFile(*t,GetValue(*t,LGAage),interval);
      writeFip(*t,interval);
    }
  }
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

//=========================================================================================================================
#define HIT_THE_FOLIAGE 1
#define NO_HIT 0
#define HIT_THE_WOOD -1


//Tassa valonlasenta voxelspacessa

// Calculation of radiation: voxel space & border forest with storing of
// Qin_stand into the segment (is otherwise as EvaluateRadiationForCfTreeSegment_1)

//This functor EvaluateRadiationForCfTreeSegment evaluates shading
//caused by all other segments on this conifer segment. This functor
//uses functor ShadingEffectOfCfTreeSegment<TS,BUD> to go through all
//segments to check the shading.

//If the attributes voxel_space and border_forest are set (then
//voxel_space != NULL), the attenuation of the bean in the voxel_space
//and border_forest is taken into consideration.

/* template <class TS, class BUD> */
/* class EvaluateRadiationForCfTreeSegment_3 { */
/* public: */
/*   EvaluateRadiationForCfTreeSegment_3(const ParametricCurve& k) : K(k), */
/*     evaluate_border(false) {} */
/*     EvaluateRadiationForCfTreeSegment_3(const ParametricCurve& k, */
/* 					VoxelSpace* vs,BorderForest* bf, bool border, */
/* 					LGMdouble a, LGMdouble b, bool sd, LGMdouble kbc, */
/* 					bool pws): */
/*       K(k),voxel_space(vs), border_forest(bf), evaluate_border(border), */
/* 	par_a(a), par_b(b), dump_self(sd), k_border_conifer(kbc), */
/* 	pairwise_self(pws){} */

/*   TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const; */
/* private: */
/*   const ParametricCurve& K; */
/*   VoxelSpace* voxel_space; */
/*   BorderForest* border_forest; */
/*   bool evaluate_border; */
/*   LGMdouble par_a, par_b; */
/*   bool dump_self; */
/*   LGMdouble k_border_conifer; */
/*   bool pairwise_self; */
/* }; */


/* This functor EvaluateRadiationForCfTreeSegment evaluates shading */
/* caused by all other segments on this conifer segment. This functor */
/* uses functor ShadingEffectOfCfTreeSegment<TS,BUD> to go through all */
/* segments to check the shading. */


/* template <class TS, class BUD> */
/*   TreeCompartment<TS,BUD>* EvaluateRadiationForCfTreeSegment_3<TS,BUD>::operator() (TreeCompartment<TS, BUD>* tc)const */
/* { */
/*   if (TS* ts = dynamic_cast<TS*>(tc)){ */
/*     SetValue(*ts, LGAQin, 0.0); */
/*     SetValue(*ts, LGAQabs, 0.0); */
/*     Radiation  conditions are not  evaluated if  the segment  has no */
/*     foliage (in practice  there would be division by  0 in computing */
/*     absorbed radiation) */
/*     if (GetValue(*ts, LGAWf) < R_EPSILON){ */
/* 	return tc; */
/*     } */

/*     Tree<TS,BUD>& tt = GetTree(*ts); */
/*     FirmamentWithMask& firmament = GetFirmament(tt); */
/*     int number_of_sectors = firmament.numberOfRegions(); */
/*     double a_dot_b = 0.0; */
/*     vector<double> radiation_direction(3); */
/*     Point middle = GetMidPoint(*ts); */

/*     vector<double> v(number_of_sectors,0.0);  */
/*     ShadingEffectOfCfTreeSegment_1<TS,BUD> s_e(ts,K,v); */
/*     This  goes  through  the  tree  and computes  shading  based  on */
/*     1)distance  light beam traverses  in foliage,  2)foliage density */
/*     and 3) inclination light beam hits the segment. */

/*     if(pairwise_self) */
/*       ForEach(tt,s_e); */
    
/*     implement  "Ip  =  Iope^(-Vp)",  s[i] =  radiation  coming  from */
/*     direction i after this */
/*     vector<double>& s = s_e.getS(); */
/*     vector<double> s(number_of_sectors, 0.0); */
/*     vector<double> qis(number_of_sectors, 0.0); */

/*     AccumulateOpticalThickness AOT(voxel_space->getXSideLength(), par_a, par_b); */

/*     for (int i = 0; i < number_of_sectors; i++){ */
/*       MJ Iop = firmament.diffuseRegionRadiationSum(i,radiation_direction); */
      

/* 	first attenuation in the voxel space */
/* 	LGMdouble transmission_voxel = 1.0; */
/* 	vector<VoxelMovement> vm; */
/* 	PositionVector dir(radiation_direction); */
/* 	voxel_space->getRoute(vm, middle, dir, K, false);  this shoud return only the "box route" */
/* 	with traveled lengths */
/* 	calculate the extinction coeffient */
/* 	LGMdouble optical_thickness = accumulate(vm.begin(),vm.end(),0.0,AOT); */
/* 	cout << "OD " << optical_thickness << endl; */
/* 	exit(0); */
/* 	cout << "Optical thickness " << optical_thickness << endl; */
/*       Vahenna oma (segmentti) vaikutus pois jos puu mukana laskuissa */
/* 	if(dump_self) { */
/* 	  LGMdouble k; */
/* 	  if(vm[0].n_segs_real > 0.0) */
/* 	    k = max(0.0,-0.014+1.056*vm[0].STAR_mean); */
/* 	  else */
/* 	    k = 0.0; */
/* 	  optical_thickness -= k * GetValue(*ts,LGAAf) * vm[0].l / voxel_space->getBoxVolume(); */

/* 	  if(optical_thickness < 0.0) */
/* 	    optical_thickness = 0.0; */
/* 	} */
/* 	if(optical_thickness > R_EPSILON){ */
/* 	  if(optical_thickness < 20.0){ */
/* 	    transmission_voxel = exp(-optical_thickness); */
/* 	  } */
/* 	  else{ */
/* 	    transmission_voxel = 0.0; */
/* 	  } */
/* 	} */
/* 	Iop *= transmission_voxel; */

/* 	then attenuation in the BorderForest */
/* 	if(evaluate_border) */
/* 	   Iop *= border_forest->getBorderForestExtinction(middle, dir, k_border_conifer); */
      
/*       qis[i] = Iop; */

/*       if(pairwise_self) { */
/* 	if (s[i] == HIT_THE_WOOD){ */
/* 	  s[i] = 0.0; */
/* 	} */
/* 	else */
/* 	  s[i] = Iop*exp(-s[i]); */
/*       } */
/*       else { */
/* 	s[i] = Iop; */
/*       }	 */

/*     } End of no_sectors ... */


/*    Total incoming radiation and radiation after stand */
/*     LGMdouble Qin_stand = accumulate(qis.begin(),qis.end(),0.0); */
/*     ts->setQinStand(Qin_stand); */

/*     MJ Q_in = accumulate(s.begin(),s.end(),0.0); */
    
/*     s contains now incoming radiation from each sector. Evaluate how */
/*     much segment absorbs from incoming radation. */
/*     LGMdouble Lk, inclination, Rfk, Ack, extinction, sfk, Ask, Wfk; */
/*     Lk = Rfk = Ack =  extinction = sfk = Ask = Wfk = 0.0; */
/*     Lk = GetValue(*ts, LGAL);   length is > 0.0, otherwise we would not bee here */
/*     Rfk = GetValue(*ts, LGARf);  Radius to foliage limit  */
/*     Wfk = GetValue(*ts, LGAWf); Foliage mass */
/*     sfk  = GetValue(tt, LGPsf); //Foliage m2/kg from tree */
/*     sfk  = GetValue(*ts, LGAsf); Foliage m2/kg from segment!!! */

/*     for (int i = 0; i < number_of_sectors; i++){ */
/*       firmament.diffuseRegionRadiationSum(i,radiation_direction); */
/*       a_dot_b = Dot(GetDirection(*ts), PositionVector(radiation_direction)); */
/*       inclination = PI_DIV_2 - acos(fabs(a_dot_b)); */

/*       Ack = 2.0*Lk*Rfk*cos(inclination) + PI_VALUE*pow(Rfk,2.0)*sin(inclination); */
/*       extinction = (double)K(inclination); */

/*       if (Ack == 0.0){ */
/* 	cout << "ERROR EvaluateRadiationForCfTreeSegment: Ack == 0 (division by 0)" */
/* 	     << endl; */
/*       } */

/*       implement I(k)p = Ip*Ask, Note  Ack must be greater than 0 (it */
/*       should if there is any foliage) */
/*       Ask = (1.0 - exp(-extinction*((sfk*Wfk)/Ack)))*Ack; */
/*       s[i] *= Ask; */
/*     } */
/*     MJ Q_abs = accumulate(s.begin(),s.end(),0.0); */
/*     SetValue(*ts, LGAQabs, Q_abs); */
/*     SetValue(*ts, LGAQin, Q_in); */
/*   } */
/*   return tc; */
/* } */


#undef HIT_THE_FOLIAGE
#undef NO_HIT
#undef HIT_THE_WOOD


#endif
