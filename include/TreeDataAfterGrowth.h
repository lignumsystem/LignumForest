#ifndef TREEDATA_AFTER_GROWTH_H
#define TREEDATA_AFTER_GROWTH_H
#include <array>
#include <vector>
#include <utility>
#include <map>
#include <SomeFunctors.h>
///\file TreeDataAfterGrowth.h
///\brief Define data structures to be able to write HDF5 data file after simulation.

///Group name for parameters
const string PGROUP("/Parameters/");
///Group name for functions in Lignum Tree class
const string TFGROUP("/TreeFunctions/");
///Group name for all functions used, i.e. all `*.fun` files found
const string AFGROUP("/AllFunctions/");
///Group name for command line values
const string CLGROUP("/CommandLine/");
///Group name for trees as xml std::string
const string TXMLGROUP("/TreeXML/");
/// HDF5 Dataset name for command line
const string COMMAND_LINE_DATASET_NAME("CommandLine");
/// HDF5 Dataset name for stand level data
const string STAND_DATA_DATASET_NAME("StandData");
/// HDF5 dataset name for center stand level data
const string CENTER_STAND_DATA_DATASET_NAME("CenterStandData");
/// HDF5 attribute name for stand data column names.
const string STAND_DATA_COLUMN_ATTRIBUTE_NAME("StandDataColumnNames");
/// Column names for stand data
const array<string,19> sdcn={
  "Year","StandArea","N_trees","10000*N_trees/StandArea","Dbase_mean","Dbase_min","Dbase_max","Dbh_mean","Dbh_min","Dbh_max",
  "H_mean","H_min","H_max","StandBasalArea","StandBasalAreaCrownBase","StandStemVol","LAI","Stand_Wf",
  "CrownLimit_mean"
};
/// It should be possible to use {} initialization throughout C++
/// Using vectors in Lignum implementation
const vector<string> STAND_DATA_COLUMN_NAMES(sdcn.begin(),sdcn.end());
/// HDF5 Dataset name Tree level data
const string TREE_DATA_DATASET_NAME("ForestTreeData");
/// HDF5 attribute name for tree level data (i.e. the dataset) column names 
const string TREE_DATA_COLUMN_ATTRIBUTE_NAME("TreeDataColumnNames");
/// Column names for the HDF5 file for tree data.
/// \sa GrowthLoop::collectDataAfterGrowth TREE_DATA_COLUMN_NAMES 
const array<string,51> tdcn={
  "TreeId","X","Y","Z","TreeNseg","TreeCrownVol","TreeH","TreeDBase","TreeDbh","TreeDCrownBase",
  "TreeHCrownBase","TreeAsBase","TreeAsDbh", "TreeAsCrownBase","TreeAf","AxisVol","TreeP","TreeM","Mr_prev","M_above",
  "Ms","Mf","Wf","Wf_new","Ws","Ws_old","Ws_D_growth","Ws_new","Ws_D_growth+Ws_new", "Wr","Wr_new",
  "QinTop","QinMax","QinTop/QinMax","Qabs","Qabs/(DiffBallSensor*TreeAf)","Wf_P","TreeP/Wf_P","ASeg0","W","Wstem","Wbranch","Ws_stem","Nsegment",
  "MeanBranch_SumD^2","MeanBranch_SumL","MeanBranch_SumD^2*L","MeanBranch_SumD^2*L/SumD^2","MeanBranch_Nbranch","MeanBranch_SumL/Nbranch","lambda"
};
/// Column names for the HDF5 file tree data.
/// -# TreeId:   Unique number of the tree in the Forest
/// -# X:   The x,y,z coordinates define the location of the tree
/// -# Y:   y xoordinate
/// -#  Z:   z coordinate
/// -#  TreeNseg:   Number of segments in a tree
/// -#  TreeCrownVol:   Crown volume \sa CrownVolume
/// -#  TreeH:   Tree height
/// -#  TreeDBase:   Tree base diameter
/// -#  TreeDbh:   Tree diameter breast height (D 1.3m)
/// -#  TreeDCrownBase:   Tree diameter crown base
/// -#  TreeHCrownBase:   Tree crown base height
/// -#  TreeAsBase:   Sapwood area at base
/// -#  TreeAsDbh:   Sapwood area Dbh (D1.3 m)
/// -#  TreeAsCrownBase:   Sappwood area at crown base
/// -#  TreeAf:   Foliage area in a tree
/// -#  AxisVol:   Main axis volume \sa MainAxisVolume
/// -#  TreeP:   Tree photosynthesis (*before* new growth tree segments)
/// -#  TreeM:   Tree respiration (*before* new growth)
/// -#  Mr_prev:   Root respiration *before* new growth
/// -#  M_above:   Above ground respiration:   TreeM-Mr_prev
/// -#  Ms:   Sapwood respiration (*before* new growth)
/// -#  Mf:   Folaige respiration (*before* new growth)
/// -#  Wf:   Tree foliage mass
/// -#  Wf_new:   Foliage mass in *new* segments
/// -#  Ws:   All sapwood mass (Ws_new+Ws_old)
/// -#  Ws_old:   Sapwood mass in *old* segments after growth
/// -#  Ws_D_growth:   Sapwood mass for new diameter growth:   Ws_old - sapwood_after_senescense
/// -#  Ws_new:   Sapwood in *new* segments.
/// -#  Ws_D_growth+Ws_new:   Sapwood in growth
/// -#  Wr:   Root mass
/// -#  Wr_new:   New root mass (required by foliage)
/// -#  QinTop:   Incoming radiation at the top of the tree
/// -#  QinMax:   Max Qin in the forest
/// -#  QinTop/QinMax:   Relative Qin for the tree 
/// -#  Qabs:   Absorbed radiation in a tree
/// -#  Qabs/(DiffBallSensor*TreeAf):   Radiation efficiency
/// -#  Wf_P:   foliage that photosynthesised
/// -#  TreeP/Wf_P:   P/Wf ratio
/// -#  ASeg0:   Surface area of new segments
/// -#  W:   Tree Wood mass (sapwood+heartwood)
/// -#  Wstem:   Wood mass main axis
/// -#  Wbranch:   Wood mass in branches (W-Wstem)
/// -#  Ws_stem:   Sapwood in main axis
/// -#  Nsegment:   Number of tree segments in a tree
/// -#  MeanBranch_SumD^2:   Sum of branch diameters squared (main axis branches only)
/// -#  MeanBranch_SumL:   Sum of branch lengths (branch length:   sum of segment lengths in the main axis of a branch)
/// -#  MeanBranch_SumD^2*L:   Sum of each branch D^2*L
/// -#  MeanBranch_SumD^2*L/SumD^2:   Mean branch length (Branch lengths weighted by D^2)
/// -#  MeanBranch_Nbranch:   Sum of branches in the main axis.
/// -#  MeanBranch_SumL/Nbranch:   Mean branch length (SumL/Nbranch)
/// -#  lambda:   Lambda s.t. P-M=G(lambda)
/// \note Technically C++ standard defines vector initialization as with array.
/// It seems not all compilers have implemented it yet. \sa tdcn.
const vector<string> TREE_DATA_COLUMN_NAMES(tdcn.begin(),tdcn.end());

///Tree parameters for the HDF5 file
///\sa LGMPD
const array<string,19> tree_param_names={
  "LGPaf","LGPapical","LGPar","LGPlen_random","LGPLmin","LGPlr","LGPmf","LGPmr","LGPms",
  "LGPna","LGPnl","LGPpr","LGPq", "LGPrhoW","LGPsf","LGPsr","LGPss","LGPxi","LGPzbrentEpsilon"
};
///It should be possible to use {} initialization throughout C++, compilers need to be implemented.
///Using vectors in Lignum implelementations
const vector<string> TREE_PARAMETER_NAMES(tree_param_names.begin(),tree_param_names.end());
///Dataset name for tree parameters
const string TREE_PARAMETER_DATASET_NAME("TreeParameters");
///Attribute name for HDF5 parameter dataset, i.e. column names for data frame  
const string TREE_PARAMETER_ATTRIBUTE_NAME("TreeParameterNames");
///Column names for functions used in simulations
const array<string,2> fn_columns={"X","F(X)"};
const vector<string> TREE_FN_COLUMN_NAMES(fn_columns.begin(),fn_columns.end());
///Attribute name for function datasetes, i.e. column names for data frames 
const string TREE_FN_ATTRIBUTE_NAME("ColumnNames");
///Vector of functions to write to HDF5 file
const array<LGMF,7> fna={LGMAL,LGMFM,LGMIP,LGMLONB,LGMNB,LGMVI,LGMVIONB};
const vector<LGMF> FN_V(fna.begin(),fna.end());
const array<string,7> fna_str={"LGMAL","LGMFM","LGMIP","LGMONB","LGMNB","LGMVI","LGMVIONB"};
const vector<string> FNA_STR(fna_str.begin(),fna_str.end());

#endif
