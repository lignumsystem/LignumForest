//Include Lignum implementation 
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <Lignum.h>
#include <Bisection.h>
#include <Shading.h>
#include <TreeLocations.h>
#include <CopyDumpCfTree.h>
#include <ForestFunctor.h>
//Include the implementation of the tree segment and bud
#include <ScotsPine.h>

#include <VoxelSpace.h>



//Impelements VisualizeLGMTree
/*
#include <GLSettings.h>
#include <OpenGLUnix.h>
#include <LGMVisualization.h>
*/

//Includes all kinds of stuff, turtle graphics etc.
#include <lengine.h>

//and for pine, see also pine9bp.L in lsys.
namespace Pine{
#include <LSystem.h>

}

#include <HarvestForestStand.h> 
#include <SomeFunctors.h>


int main(int argc, char** argv)
{
  int iterations = 0;//Number of iterations
  double gap_radius = 0.0; //Gap surrounding the tree
  double x_coord = 0.0;//Locations of the tree
  double y_coord = 0.0;
  string metafile;//File denoting parameters and functions
  string voxelfile;//Voxel space
  string prefix;//File prefix for output files, e.g. at CSC /wrk/jperttun
  if (argc < 7){
    cout << "Usage: ./pine Iterations GapRadius X Y MetaFile.txt ";
    cout << "VoxelSpace.txt [-startVoxCalc <ii> -viz -toFile <filename>] -treeDist <dist> -numParts <parts> ";
    cout << "-prefix <directory_name>" << endl;
    return 0;
  }
  else{
    iterations = atoi(argv[1]);
    gap_radius = atof(argv[2]);
    x_coord =  atof(argv[3]);
    y_coord =  atof(argv[4]);
    metafile = string(argv[5]);
    voxelfile = string(argv[6]);
  }

  //Check if output to file is required. Open the file here and
  //truncate, hence results of one run apper in the file (that can be
  //used by gnuplot quick-plot.gpl

  string clarg;
  bool toFile = ParseCommandLine(argc,argv,"-toFile", clarg);
  ofstream ff(clarg.c_str() , ofstream::trunc);
  if(!ff)
    toFile = true;   //If opening the output file failed
  
  int startVoxelCalculation = 10;
  if(ParseCommandLine(argc,argv,"-startVoxCalc", clarg))
    if(clarg.length() > 0)
      startVoxelCalculation = atoi(clarg.c_str());

  double tree_distance = 1.0; //Minimum distance between trees in meters
  if (ParseCommandLine(argc,argv,"-treeDist", clarg))
    if(clarg.length() > 0)
      tree_distance = atof(clarg.c_str());

  if (ParseCommandLine(argc,argv,"-prefix", clarg))
    if(clarg.length() > 0)
      prefix = clarg;

  int num_parts = 1; //Number  parts a  segment is  divided  into when
		     //assigning foliage and computing Qabs.
  if (ParseCommandLine(argc,argv,"-numParts", clarg))
    if(clarg.length() > 0)
      num_parts = atoi(clarg.c_str());

  int init = -1;//random number generator 
  ran3(&init);

  //Create the forest: trees and their l-systems  
  vector<Pine::LSystem<ScotsPineSegment,
    ScotsPineBud,PBNAME,PineBudData>*> plv;
  vector<Tree<ScotsPineSegment,
    ScotsPineBud>*> pinev;
						   
  InitializeTree<ScotsPineSegment,ScotsPineBud> init_pine(metafile);

  //Read voxel space dimension
  ifstream vf(voxelfile.c_str());
  int vx,vy,vz; vx = vy = vz = 0;
  LGMdouble s1,s2,s3; s1 = s2 = s3 = 1.0;
  vf >> vx >> vy >> vz >> s1 >> s2 >> s3;
  cout << "Voxel Space: " << vx << " " << vy << " " << vz << " " 
       << s1 << " " << s2 << " " << s3 << endl;

  //Generate tree positions, number of trees/ha in "stemsha.fun"
  ParametricCurve stems_ha("stemsha.fun");
  int nTrees = static_cast<int>(stems_ha(0)*(vx*vy/10000.0));
  ForestGap gap(pair<double,double>(vx/2,vy/2),gap_radius);
  vector<pair<double,double> > locations;
  //leave 3m distance to the forest edge (i.e.  ensure the tree crowns
  //are i n  the forest)
  double edge_width = 3.0;
  GenerateLocations(nTrees,edge_width,edge_width,vx-edge_width,
		    vy-edge_width,tree_distance,gap,locations);
  cout << "Number of locations: " << nTrees << " " << locations.size() 
       << " Density/ha wanted: " << stems_ha(0)
       << " Density/ha created: " << locations.size()*10000.0/(vx*vy) << endl; 
  cout << " Gap radius: "  << gap_radius << " Minimum tree distance: " << tree_distance  
       << " Segment parts " << num_parts << endl; 
  cout << " Tree Location: " << x_coord << " " << y_coord << endl;

  ostringstream locations_file;
  locations_file << prefix << "Locations.txt";
  ofstream fs(locations_file.str().c_str());
  //How to do this with STL copy algorithm with 'cout'??
  for(int i = 0; i < locations.size(); i++){
    fs << locations[i].first << " " << locations[i].second << endl;
  }

   //Create the forest
  for (int i = 0; i < locations.size(); i++){
    Point p(locations[i].first,locations[i].second,0);
    PositionVector up(0,0,1);
    Tree<ScotsPineSegment,ScotsPineBud> *t = new Tree<ScotsPineSegment,ScotsPineBud>(p,up);
    init_pine.initialize(*t);
    pinev.push_back(t);
    Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData>* l = 
      new Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData>();
    plv.push_back(l);
  }
  //Create the tree to be studied  in the middle of the forest, put it
  //last in the vector
  Tree<ScotsPineSegment,ScotsPineBud> *pine1 = 
    new Tree<ScotsPineSegment,ScotsPineBud>(Point(vx/2.0,vy/2.0,0),PositionVector(0,0,1));
  init_pine.initialize(*pine1);
  pinev.push_back(pine1);
  //Do the same for its L-system
  Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData>* l = 
      new Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData>();
  plv.push_back(l);
  
  //Now create the voxel space (we have the tree to get the Firmament)
  //vx,vy,vz define the voxel space  dimensions, s1, s2, s3 define the
  //voxel box dimensions
  VoxelSpace vs(Point(0,0,0),Point(vx,vy,vz),
                s1,s2,s3,
                static_cast<int>(vx/s1),static_cast<int>(vy/s2),static_cast<int>(vz/s3),
                GetFirmament(*pine1));

  //Expand axioms, first initial structure
  for (int i = 0; i < plv.size(); i++){
    plv[i]->start();
    plv[i]->lstringToLignum(*pinev[i],1,PBDATA);
    //The initial foliage and based on it the root mass
    double wf = 0.0;
    wf = Accumulate(*pinev[i],wf,CollectNewFoliage<ScotsPineSegment,ScotsPineBud>());
    //Initial root mass
    SetValue(*pinev[i],TreeWr,GetValue(*pinev[i],LGPar)*wf);
  }
  cout << "Init done" << endl;  

  //The growth loop
  for (int iter=0; iter < iterations; iter++)
  {
    cout << "Iter: " << iter << endl;
    //Resetting light calculations and needle masses, needle areas and stars. 
    //If I forget to reset something, please add!!
    vs.reset();
    //Single  tree simulation:  find the  bounding box  of  the single
    //tree,  move  the  voxelspace  so  that  its  lower  left  corner
    //coincides with the lower left corner of the bounding box, resize
    //the  voxel space  to have  the size  of the  bounding  box. (See
    //SugarMaple)

    //The star  mean currently as  a function of segment  age. Another
    //possibility is a function of its inclination.
    SetStarMean<ScotsPineSegment,ScotsPineBud> setstar(ParametricCurve("fstarmean.fun"));
    for_each(pinev.begin(),pinev.end(),ForestSetStar(setstar));
    //Harvesting the stand, stems_ha as  a function of d1.3, then dump
    //the trees into the voxel space
    double d13cm = GetValue(*pine1,LGADbh)*100.0;//meters to cm
    double trees_left = stems_ha(d13cm)*vx*vy/10000.0;
    //We   should   be   careful   in  the   "stemsha.fun"   so   that
    //trees_left/locations.size() never becomes greater than 1
    double p_remove = max(1.0 - (trees_left/locations.size()),0.0);
    HarvestForestStand(locations,pinev,plv,p_remove);
    //Now  dump the trees
    for_each(pinev.begin(),pinev.end(),ForestDumpFoliage(vs,num_parts));
    cout << "Foliage dump done" <<endl; 
    RemoveDeadTrees(locations,pinev,plv);
    cout << "Removed dead trees done" <<endl;
    //Calculation of radiation conditions, either pairwise comparison or
    //voxelspace calculations
    //pairwise  comparison
    if (iter < startVoxelCalculation){
      ParametricCurve K(0.2);
      ForestDiffuseRadiation fdr(K);
      //fdr uses DiffuseForestRadiation<ScotsPineSegment,ScotsPineBud>
      //FRad(Hc,Af,K).
      for_each(pinev.begin(),pinev.end(),fdr);
      cout << "single trees " << endl;
    }
    //use voxel space
    else{
      //true means the box shades itself
      vs.calculateTurbidLight(true);
      cout << "calculate Turbid light done" << endl;
      for_each(pinev.begin(),pinev.end(),SetForestTreeQabs(vs,num_parts));  
      cout << "set forest qabs done" << endl;
    }
    cout << "Calculated light"  <<endl;
    //Photosynthesis
    for_each(pinev.begin(),pinev.end(),ForestPhotosynthesis());
    //Respiration includes foliage, sapwood and roots
    for_each(pinev.begin(),pinev.end(),ForestRespiration());
    cout << "P and M done"  <<endl;
    //To print consistent production  and respirations we must collect
    //masses now before senescense and new growth
    double wroot = GetValue(*pine1,TreeWr);  //root mass
    //Collect foliage 
    LGMdouble wfoliage = 0.0;
    Accumulate(*pine1,wfoliage,CollectFoliageMass<ScotsPineSegment,ScotsPineBud>());
    //Collect sapwood mass
    LGMdouble wsapwood = 0.0;
    Accumulate(*pine1,wsapwood,CollectSapwoodMass<ScotsPineSegment,ScotsPineBud>());
    
    //After respiration the senescence, includes root senescence
    for_each(pinev.begin(),pinev.end(),ForestAging());
    
    //This first derive() creates the new segments, whose lengths will be iterated
    cout << "Begin new segments" <<endl;
    for (int i = 0; i < plv.size(); i++){
      plv[i]->derive();
      //This endEach will reset the mode back to 0
      plv[i]->endEach();
      plv[i]->lstringToLignum(*(pinev[i]),1,PBDATA);
    }
    //This endEach will set the mode  to 1
    plv[0]->endEach();
    cout << "End new segments" <<endl;
    //Pass the qin to newly created segments
    for_each(pinev.begin(),pinev.end(),ForestForwardQin());
    //The vigour index
    for_each(pinev.begin(),pinev.end(),ForestVigourIndex());
    
    if (iter < startVoxelCalculation){
      //Ball sensor 
      for_each(pinev.begin(),pinev.end(),ForestPairwiseQinMax());
    }
    else{
      //In the voxel space, at the moment the tree itself
      for_each(pinev.begin(),pinev.end(),ForestVoxelSpaceQinMax());
    }
    cout << "Begin ForestPathLength" <<endl;
    for_each(pinev.begin(),pinev.end(),ForestPathLength());
    cout << "End ForestPathLength" <<endl;
    //ForestAllocateGrowth  does: Allocation, KillBudsAfterAllocation,
    //TreeDiameterGrowth and updates root mass
    cout << "Begin carbon allocation" <<endl;
    for_each(pinev.begin(),pinev.end(),ForestAllocateGrowth());
    cout << "End carbon allocation" <<endl;
    RemoveDeadTrees(locations,pinev,plv);

    //Update L-string, pass the state of the bud to control branching
    //Before derive  pass the  foliage mass of  the mother  segment to
    //terminating buds
    for_each(pinev.begin(),pinev.end(),ForwardForestWf());
    cout << "Begin L-string update" <<endl;
    for (int i = 0; i < plv.size(); i++){
      plv[i]->lignumToLstring(*pinev[i],1,PBDATA);
    }
    cout << "End L-string update" <<endl;
    //create new buds as the function of the mother segment foliage mass
    cout << "Begin new buds" <<endl;
    for (int i = 0; i < plv.size(); i++){
      plv[i]->derive();
      //The endEach will reset the mode to 1
      plv[i]->endEach();
      plv[i]->lstringToLignum(*pinev[i],1,PBDATA);
    }
    //The endEach will set the mode to 0
    plv[0]->endEach();
    cout << "End new buds" << endl;
    Axis<ScotsPineSegment,ScotsPineBud>& axis = GetAxis(*pine1);
    //Collect foliage after growth
    LGMdouble wfaftergrowth = 0.0;
    Accumulate(*pine1,wfaftergrowth,CollectFoliageMass<ScotsPineSegment,ScotsPineBud>());
    //Collect sapwood mass after growth
    LGMdouble wsaftergrowth = 0.0;
    Accumulate(*pine1,wsaftergrowth,CollectSapwoodMass<ScotsPineSegment,ScotsPineBud>());
    list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>& ls = GetTreeCompartmentList(axis);
    //Mean branch length
    summing bs;
    bs.d2 = bs.d2l = bs.lsum = 0.0;
    bs.n_br = 0;
    bs = Accumulate(*pine1, bs, Branchmeans() );

    //The Qin at the top
    LGMdouble qintop1 = 0.0;
    GetTopQin<ScotsPineSegment,ScotsPineBud> getTopQin1;
    qintop1 = accumulate(ls.begin(),ls.end(),qintop1, getTopQin1);

    //Diameter and heigth at the crown base.
    DCLData dcl;
    AccumulateDown(*pine1,dcl,AddBranchWf(),DiameterCrownBase<ScotsPineSegment,ScotsPineBud>());

    //Crown volume
    CrownVolume<ScotsPineSegment,ScotsPineBud> cv(0.30);
    double cvol = cv(*pine1);

    int nsegment = 0;
    nsegment = AccumulateDown(*pine1,nsegment,
			      CountTreeSegments<ScotsPineSegment,ScotsPineBud>()); 
    //There are  some gnuplot  scripts that can  plot data  files from
    //this output.   Please just add  additional output to the  end if
    //needed (please do not  'break' those scripts).  Please note: Due
    //to logic  in the  main loop the  production and  respiration are
    //done before  new growth and  senescence. The masses  plotted are
    //after senescence and new growth.
    cout << "To file" << endl;
    if(toFile) {
      double qbs = 0.0;
      double qabs= Accumulate(*pine1,qbs,CollectQabs<ScotsPineSegment,
			    ScotsPineBud>());
      double tla0 = 0.0;
      double treeAf = Accumulate(*pine1,tla0,CollectFoliageArea<
				 ScotsPineSegment,ScotsPineBud>());

      ff <</*1 */ left << setw(6) << setfill(' ') << iter+1 << " "
	 <</*2 */ left << setw(8) << locations.size()*10000.0/(vx*vy) << " " 
	 <</*3 */ left << setw(8) << GetValue(*pine1,LGAH) << " " 
	 <</*4 */ setw(9) << qintop1 << " " 
	 <</*5 */ setw(12) << qintop1/GetValue(*pine1,TreeQinMax) << " " 
	 <</*6 */ setw(11) << GetValue(*pine1,LGADbase) << " "
	 <</*7 */ setw(11) << GetValue(*pine1,LGADbh) << " "
	 <</*8 */ setw(11) << dcl.DCrownBase() << " " 
	 <</*9 */ setw(11) << dcl.HCrownBase() << " " 
	 <</*10*/ setw(11) << GetValue(*pine1,TreeP) << " " //Production before new growth
	 <</*11*/ setw(11) << GetValue(*pine1,TreeM) << " "//Respiration before new growth: foliage+sapwood+roots
	 <</*12*/ setw(12) << wfaftergrowth << " " //Foliage mass after new growth
	 <</*13*/ setw(14) << bs.d2l/bs.d2 << " "  //Branch means
	 <</*14*/ setw(14) << bs.lsum/(double)bs.n_br  << " "
	 <</*15*/ setw(11) << GetValue(*pine1,TreeWr) << " "     //Root mass after new growth
	 <</*16*/ setw(11) << GetValue(*pine1,LGPmr)*wroot << " "//Root respiration before new growth
	 <</*17*/ setw(11) << GetValue(*pine1,TreeM) - GetValue(*pine1,LGPmr)*wroot << " "//Above ground respiration 
	 <</*18*/ setw(11) << wsaftergrowth  << " "           //Sapwood mass after new growth
	 <</*19*/ setw(11) << GetValue(*pine1,LGPms)*wsapwood << " " //Sapwood respiration before new growth
	 <</*20*/ setw(11) << GetValue(*pine1,LGPmf)*wfoliage << " " //Foliage respiration before new growth
	 <</*21*/ setw(11) << cvol << " "  // Crown volume
	 <</*22*/ setw(11) << 0.0 << " "//Sum[(LeafA+NeedleA)/Vbox]/Nboxes
	 <</*23*/ setw(11) << nsegment << " " //Number of segments
	 <</*24*/ setw(11) << qabs << " "     // Qabs  
	 <</*25*/ setw(11) << qabs/(GetFirmament(*pine1).
				    diffuseBallSensor()*treeAf) << " "  //Rad eff.
	 <</*26*/ setw(11) << wfoliage << " " //Foliage that photosynthesized
	 <</*27*/ setw(11) << GetValue(*pine1,TreeP)/wfoliage << " "  //The ubiquitous P/Wf
	 << endl; //Qabs
    }
    ostringstream summary_file_name;
    summary_file_name << prefix << "Trees-" << iter+1 << ".txt";
    ofstream summary_file(summary_file_name.str().c_str());
    cout << "To summary file " << summary_file_name.str() << endl;
    for_each(pinev.begin(),pinev.end(),ForestPrintSummary(summary_file));
    
    cout << "To crown limit file" <<endl;
    //Collect and print the two lists of foliage masses and their heights 
    //in the main axis.  Print and assess the crown base afterwards.
    CrownLimitData cld;
    AccumulateDown(*pine1,cld,AddCrownLimitData(),
		   CollectCrownLimitData<ScotsPineSegment,ScotsPineBud>());

    ostringstream crown_limit_file;
    crown_limit_file << prefix << "CrownLimit-" << x_coord << "-" << y_coord << "-" << iter+1 << ".txt";
    //File is CrownLimit+iter+.txt, e.g., "CrownLimit-25.5-22.0-10.txt"
    ofstream cl_file(crown_limit_file.str().c_str());
    const list<pair<double,double> >& hwf_ls = cld.WfHList();
    const list<pair<double,double> >& dh_dwf_ls = cld.dHdWfList();
    const list<pair<double,double> >& h_qabs_ls = cld.HQabsList();
    const list<pair<double,double> >& dh_dqabs_ls = cld.dHdQabsList();
    //How to do this with copy algorithm as with 'cout'??
    list<pair<double,double> >::const_iterator hwf_it = hwf_ls.begin();
    list<pair<double,double> >::const_iterator dh_dwf_it = dh_dwf_ls.begin();
    list<pair<double,double> >::const_iterator h_qabs_it =  h_qabs_ls.begin();
    list<pair<double,double> >::const_iterator dh_dqabs_it = dh_dqabs_ls.begin();
    cl_file << setfill(' ') 
	    << setw(11) << "H" << " " <<  setw(11) << "Wf" << " "
	    << setw(11) << "dH" << " " <<  setw(11) << "dWf/dH" << " " 
	    << setw(11) << "H" << " " <<  setw(11) << "Qabs" << " " 
	    << setw(11) << "dH" << " " <<  setw(11) << "dQabs/dH" << endl; 
    while (hwf_it != hwf_ls.end()){
      cl_file << setfill(' ')
	      << setw(11) << (*hwf_it).first << " " 
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
    }
    ostringstream voxel_space_file("VoxelSpace-");
    voxel_space_file << "VoxelSpace-" << x_coord << "-" << y_coord << "-" << iter+1 << ".txt";
    //File is VoxelSpace+location+iter.txt, e.g., "VoxelSpace-25.5-22.0-10.txt"
    //vs.writeVoxelBoxesToGnuPlotFile(voxel_space_file);
    cout << "Begin prune" <<endl;
    for (int i = 0; i < plv.size(); i++){
      plv[i]->prune(*pinev[i]);
    }
    cout << "End prune" <<endl;
  }   // END OF ITERATION END OF ITERATION END OF ITERATION END OF ITERATION


  //Clean up.
  cout << "Growth end" << endl;

  //Close result file if was open
  if(toFile)
    ff.close();
  //Possible clean up in L-system
  for (int i = 0; i < plv.size(); i++){
    plv[i]->end();
  }  
  //Print  Qin and  (x,y,z)  for each  segment  of age  0,1  or 2.  In
  //addition  to Qin,  this gives  you  the idea  of the  size of  the
  //branches
  ForEach(*pine1,PrintSegmentQin<ScotsPineSegment,ScotsPineBud>());
  LGMdouble wf1 = 0.0;
  cout << "Wf: " << Accumulate(*pine1,wf1,
       		       CollectFoliageMass<ScotsPineSegment,ScotsPineBud>()) << endl;
  summing bs;
  bs = Accumulate(*pine1, bs, Branchmeans());
  if(bs.d2 > 0.0)
    cout << "Averaged2 branch length, m: " << bs.d2l/bs.d2 << endl;
  else
    cout << "Somethin wrong in d2 branch averaging!" << endl;
  

  if(bs.n_br > 0)
    cout << "Average branch length, m: " << bs.lsum/(double)bs.n_br
	 << "   no. branches: " << bs.n_br << endl;
  else
    cout << "Somethin wrong in d2 branch averaging!" << endl;

  int n = 0;
  cout << "Segments: " << Accumulate(*pine1,n, 
				     CountTreeSegments<ScotsPineSegment,
				     ScotsPineBud>()) << endl;
  LGMdouble qabs = 0.0;
  cout << "Tree       Qabs: " 
       <<  Accumulate(*pine1,qabs,CollectQabs<ScotsPineSegment,ScotsPineBud>())<<endl;
// cout << "VoxelSpace Qin : " <<  vs.getQin() << endl;
  ForEach(*pine1,MoveTree<ScotsPineSegment,ScotsPineBud>(Point(-GetPoint(*pine1).getX(),
					      -GetPoint(*pine1).getY(),0)));
 
  //Print the locations. Generate the file name Locations<n>.txt where
  //n is the value of 'iter'

  locations_file.clear();
  locations_file << prefix <<"Locations-" << iterations << ".txt";
  cout << "Locations: " << locations_file.str().c_str() << endl;
  ofstream ofs(locations_file.str().c_str());
  //How to do this with copy as with 'cout'??
  for(int i = 0; i < locations.size(); i++){
    ofs << locations[i].first << " " << locations[i].second << endl;
  }

  /*   
  if(CheckCommandLine(argc,argv,"-viz")) {

    LGMVisualization viz;
    viz.InitVisualization();

  
    // textures 512x512
    viz.AddCfTree(*pine1, "Manty.bmp", "neulaset5.tga");
    //viz.makeDisplayLists();
    // viz.hello();
    float th = (float)GetValue(*pine1,LGAH);
    cout << th << endl;
    viz.ResetCameraPosition(th);
    viz.SetMode(WIREMODEL);
    viz.ResetCameraPosition(GetValue(*pine1,LGAH));
    viz.StartVisualization();
  }
  */
  return 0;
 
}
