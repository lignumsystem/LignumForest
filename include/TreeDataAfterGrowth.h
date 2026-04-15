#ifndef TREEDATA_AFTER_GROWTH_H
#define TREEDATA_AFTER_GROWTH_H
#include <array>
#include <vector>
#include <utility>
#include <map>
#include <SomeFunctors.h>
/// \file TreeDataAfterGrowth.h
/// \brief Names and data structures for HDF5 tables after simulation.
///
/// \defgroup HDF5FILES HDF5 file constants
/// \brief Constants for HDF5 file names, groups, datasets and dataset attributes
/// @{
/// \defgroup HDF5GROUPNAMES HDF5 file group names
/// \defgroup HDF5COMMANDLINEDATASET HDF5 file command line dataset
/// \defgroup HDF5STANDDATASET HDF5 file stand dataset
/// \defgroup HDF5TREEDATASET HDF5 tree dataset
/// \defgroup HDF5DEADTREEDATASET HDF5 dead tree dataset
/// \defgroup HDF5PARAMETERDATASET HDF5 file parameter dataset
/// \defgroup HDF5FUNCTIONDATASET HAD5 file function dataset
/// \defgroup HDF5VOXELSPACESIZEDATASET HDF5 file voxel space size dataset
/// \defgroup HDF5VOXELSPACEDATADATASET HDF5 file voxel space data dataset
/// @}
namespace LignumForest{
  /// \addtogroup HDF5GROUPNAMES
  ///@{
  ///\name HDF5 file group constants
  ///Names and file prefixes for HDF5 groups
  ///@{
  ///\brief Group name for parameters
  const string PGROUP("/Parameters/");
  ///\brief Group name for functions in Lignum Tree class
  const string TFGROUP("/TreeFunctions/");
  ///\brief Group name for all functions used, i.e. all `*.fun` files found
  const string AFGROUP("/AllFunctions/");
  ///\brief Group name for command line values
  const string CLGROUP("/CommandLine/");
  ///\brief Group name for trees as xml std::string
  const string TXMLGROUP("/TreeXML/");
  ///\brief Group name for Tree parameter files
  const string PFILEGROUP("/AllParameters/");
  ///\brief Group name for MetaFile files
  const string ALLMETAFILEGROUP("/AllMetaFiles/");
  ///\brief Group name for parameter files
  const string ALLPARAMFILEGROUP("/AllParameterFiles/");
  ///\brief Group name for function files
  const string ALLFNFILEGROUP("/AllFunctionFiles/");
  ///\brief Group name for Firmament
  const string FIRMAMENTGROUP("/Firmament/");
  ///\brief Group prefix for the VoxelSpaces 
  const string VOXELSPACEGROUP("/VoxelSpace/");
  ///\brief Group name for VoxelSpace sizes during sikmulation
  const string VOXELSPACESIZESGROUP("/VoxelSpaceSizes/");
  ///\brief HDF5 file name prefix for XML trees.
  const string TREEXML_PREFIX("TreesXML_");
  ///\brief HDF5 file name prefix for VoxelSpace data
  const string VOXELSPACE_PREFIX("VoxelSpaces_");
  ///@}
  ///@}
  /// \addtogroup HDF5COMMANDLINEDATASET
  ///@{
  ///\name HDF5 file command line dataset constants
  ///Command line dataset name 
  ///@{
  ///\brief HDF5 Dataset name for command line
  const string COMMAND_LINE_DATASET_NAME("CommandLine");
  ///@}
  ///@}
  /// \addtogroup HDF5STANDDATASET
  ///@{
  ///\name HDF5 file stand data and center stand constants
  ///Dataset, attribute and column names for stand and center stand datasets
  ///@{
  ///\brief HDF5 Dataset name for stand level data
  const string STAND_DATA_DATASET_NAME("StandData");
  ///\brief HDF5 dataset name for center stand level data
  const string CENTER_STAND_DATA_DATASET_NAME("CenterStandData");
  ///\brief HDF5 attribute name for stand data column names.
  const string STAND_DATA_COLUMN_ATTRIBUTE_NAME("StandDataColumnNames");
  ///\brief Column names for stand data
  const array<string,19> sdcn={
    "Year","StandArea","N_trees","10000*N_trees/StandArea","Dbase_mean","Dbase_min",
    "Dbase_max","Dbh_mean","Dbh_min","Dbh_max","H_mean","H_min","H_max","StandBasalArea",
    "StandBasalAreaCrownBase","StandStemVol","LAI","Stand_Wf","CrownLimit_mean"
  };
  ///\brief Column names for stand data
  const vector<string> STAND_DATA_COLUMN_NAMES(sdcn.begin(),sdcn.end());
  ///@}
  ///@}
  /// \addtogroup HDF5TREEDATASET
  ///@{
  ///\name HDF5 file tree dataset constants
  ///Dataset, attribute and column names for tree datasets
  ///@{
  ///\brief HDF5 Dataset name for individual tree data
  const string TREE_DATA_DATASET_NAME("ForestTreeData");
  ///\brief HDF5 attribute name for tree level data (i.e. the dataset) column names 
  const string TREE_DATA_COLUMN_ATTRIBUTE_NAME("TreeDataColumnNames");
  ///\brief Column names for the HDF5 file for tree data.
  /// \sa GrowthLoop::collectDataAfterGrowth TREE_DATA_COLUMN_NAMES 
  const array<string,52> tdcn={
    "TreeId","X","Y","Z","TreeNseg","TreeCrownVol","TreeH","TreeDBase","TreeDbh","TreeDCrownBase",
    "TreeHCrownBase","TreeAsBase","TreeAsDbh", "TreeAsCrownBase","TreeAf","AxisVol","TreeP","TreeM","Mr_prev","M_above",
    "Ms","Mf","Wf","Wf_new","Ws","Ws_old","Ws_D_growth","Ws_new","Ws_D_growth+Ws_new", "Wr","Wr_new",
    "QinTop","QinMax","QinTop/QinMax","Qabs","Qabs/(DiffBallSensor*TreeAf)","Wf_P","TreeP/Wf_P","ASeg0",
    "W","Wstem","Wbranch","Ws_stem","Nsegment","MeanBranch_SumD^2","MeanBranch_SumL","MeanBranch_SumD^2*L",
    "MeanBranch_SumD^2*L/SumD^2","MeanBranch_Nbranch","MeanBranch_SumL/Nbranch","lambda",
    "MaxBranch"
  };
  ///\brief Column names for the HDF5 file tree data:
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
  /// -#  MeanBranch_SumL:   Sum of branch lengths
  ///  -# (branch length: sum of segment lengths in the main axis of a branch)
  /// -#  MeanBranch_SumD^2*L:   Sum of each branch D^2*L
  /// -#  MeanBranch_SumD^2*L/SumD^2:   Mean branch length (Branch lengths weighted by D^2)
  /// -#  MeanBranch_Nbranch:   Sum of branches in the main axis.
  /// -#  MeanBranch_SumL/Nbranch:   Mean branch length (SumL/Nbranch)
  /// -#  lambda:   Lambda s.t. P-M=G(lambda)
  /// -#  MaxBranch: Max branch extension from the main stem 
  /// \remark  Technically C++ standard defines vector initialization as with array.
  /// It seems not all compilers have implemented it yet.
  /// \sa tdcn
  const vector<string> TREE_DATA_COLUMN_NAMES(tdcn.begin(),tdcn.end());
  ///@}
  ///@}
  /// \addtogroup HDF5DEADTREEDATASET
  ///@{
  ///\name HDF5 file dead tree dataset constants
  ///Dataset, attribute and column names for dead tree dataset
  ///@{
  ///\brief HDF5 dataset name for individual dead tree data.
  const string DEAD_TREE_DATA_DATASET_NAME="ForestDeadTreeData";
  ///\brief DF5 attribute name for individual dead tree  data (i.e. the dataset) column names
  const string DEAD_TREE_DATA_COLUMN_ATTRIBUTE_NAME("ForestDeadTreeDataColumnNames");
  ///\brief Column names for individual dead tree data. Expand as necessary.
  ///\remark Names as for TREE_DATA_COLUMN_NAMES.
  const array<string,5> dead_tdcn={"TreeId","X","Y","Z","AxisVol"};
  /// Column names for the HDF5 file dead tree data.
  const vector<string> DEAD_TREE_DATA_COLUMN_NAMES(dead_tdcn.begin(),dead_tdcn.end());
  ///@}
  ///@}
  /// \addtogroup HDF5PARAMETERDATASET
  ///@{
  ///\name HDF5 file tree parameter dataset constants
  ///Dataset, attribute and column names for  tree parameter dataset
  ///@{
  ///\brief Tree parameters for the HDF5 file
  ///\sa LGMPD
  const array<string,19> tree_param_names={
    "LGPaf","LGPapical","LGPar","LGPlen_random","LGPLmin","LGPlr","LGPmf","LGPmr","LGPms",
    "LGPna","LGPnl","LGPpr","LGPq", "LGPrhoW","LGPsf","LGPsr","LGPss","LGPxi","LGPzbrentEpsilon"
  };
  ///\brief Tree parameter names. 
  ///\remark It should be possible to use {} initialization throughout C++, compilers need to be implemented.
  const vector<string> TREE_PARAMETER_NAMES(tree_param_names.begin(),tree_param_names.end());
  ///\brief Column names for tree parameters
  ///@{
  const array<string,2> p_columns={"Name","Value"};
  const vector<string> TREE_PARAMETER_COLUMN_NAMES(p_columns.begin(),p_columns.end());
  ///@}
  ///\brief Attribute name for tree parameter attributes
  const string TREE_PARAMETER_FILE_ATTRIBUTE_NAME("ColumnNames");
  ///\brief Dataset name for tree parameters
  const string TREE_PARAMETER_DATASET_NAME("TreeParameters");
  ///\brief Attribute name for HDF5 parameter dataset, i.e. column names for data frame  
  const string TREE_PARAMETER_ATTRIBUTE_NAME("TreeParameterNames");
  ///@}
  ///@}
  /// \addtogroup HDF5FUNCTIONDATASET
  ///@{
  ///\name HDF5 file function dataset constants
  ///Dataset, attribute and column names for tree function dataset  
  ///@{
  ///\brief Column names for functions used in simulations
  const array<string,2> fn_columns={"X","F(X)"};
  const vector<string> TREE_FN_COLUMN_NAMES(fn_columns.begin(),fn_columns.end());
  ///\brief Attribute name for function datasets, i.e. column names for data frames 
  const string TREE_FN_ATTRIBUTE_NAME("ColumnNames");
  ///\brief Function names as enumerations
  const array<LGMF,7> fna={LGMAL,LGMFM,LGMIP,LGMLONB,LGMNB,LGMVI,LGMVIONB};
  ///\brief Vector name for the user
  const vector<LGMF> FN_V(fna.begin(),fna.end());
  ///\brief String names for dataset attribute
  const array<string,7> fna_str={"LGMAL","LGMFM","LGMIP","LGMONB","LGMNB","LGMVI","LGMVIONB"};
  ///\brief Vector name for the user
  const vector<string> FNA_STR(fna_str.begin(),fna_str.end());
  ///@}
  ///@}
  /// \addtogroup HDF5VOXELSPACESIZEDATASET
  ///@{
  ///\name HDF5 file voxel space sizes dataset constants
  ///Dataset, attribute and column names for VoxelSpace size development dataset.
  ///@{
  ///\brief VoxelSpace sizes dataset name
  const string VOXELSPACESIZES_DATASET_NAME("VoxelSpaceSizesData");
  ///\brief HDF5 attribute name (i.e. column names) for VoxelSpace size developnent
  const string VOXELSPACESIZES_ATTRIBUTE_NAME("VoxelSpaceSizesColumnNames");
  ///\brief LL = Lower Left corner point, UR = Upper Right corner point
  const array<string,12> vs_sizes_column_names={"Year","LLX","LLY","LLZ","URX","URY","URZ","Area","Width","Length","Height","Width_x_Length"};
  ///\brief Vector name for the user
  const vector<string> VS_SIZES_COLUMN_NAMES(vs_sizes_column_names.begin(),vs_sizes_column_names.end());
  ///@}
  ///@}
  /// \addtogroup HDF5VOXELSPACEDATADATASET
  ///@{
  ///\name HDF5 file VoxelSpace datasets constants
  ///Dataset, attribute and columns names for VoxelSpace datasets
  ///@{
  ///\brief VoxelSpace dataset name
  const string VOXELSPACE_DATA_DATASET_NAME("VoxelSpaceData");
  ///\brief VoxelSpace columns attribute name
  const string VOXELSPACE_ATTRIBUTE_NAME("VoxelSpaceColumnNames");
  ///\brief VoxelSpace data collected VBDATA dimension attribute name
  ///\sa VS_COLUMN_NAMES
  ///\sa VB_DATA_COLUMN_NAMES
  const string VOXELBOX_DATA_ATTRIBUTE_NAME("VBDATAColumnNames");
  ///\brief Column names for the VoxelSpace data
  const array<string,4> vs_column_names={"X","Y","Z","VBDATA"};
  ///\brief Vector name for the user
  const vector<string> VS_COLUMN_NAMES(vs_column_names.begin(),vs_column_names.end());
  ///\brief Data collected from a voxelspace::VoxelBox in VBDATA dimension
  const array<string,10> vb_data_column_names={"STAR", "STARsum","VAL_c","VAL_b","LGAWf","LGAAf",
					      "LGAWs","LGAAWs","Nsegments","NSeg/Nparts"};
  ///\brief Vector name for the user
  const vector<string> VB_DATA_COLUMN_NAMES(vb_data_column_names.begin(),vb_data_column_names.end());
  ///\brief Enumeration for indexing
  ///
  ///Enumeration for data matrix indexing and for HDF5 dataset creation
  ///\remark NSEGREAL: Segment position in voxelspace::VoxelSpace is checked from multiple points
  enum VB_DATA_NAMES{VB_STAR,     ///< STAR value
		     VB_STARSUM,  ///< STAR sum
		     VB_VAL_c,    ///< Conifers weighted STAR
		     VB_VAL_b,    ///< Broadleaf weighted star
		     VB_LGAWf,    ///< Foliage mass
		     VB_LGAAf,    ///< Foliage area
		     VB_LGAWs,    ///< Wood mass
		     VB_LGAAWs,   ///< Wood area
		     VB_NSEG,     ///< Number of segments
		     VB_NSEGREAL  ///< NSEG/SegmentParts 
  };
  const string VB_EDGE_SIZE_NAME("VoxelEdgeSizeXYZ");
  ///@}
  ///@}
}//End namespace LignumForest
#endif
