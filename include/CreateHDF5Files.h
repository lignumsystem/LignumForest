#ifndef CREATEHDF5FILES_H
#define CREATEHDF5FILES_H
#include <deque>
#include <LGMHDF5File.h>
#include <TreeDataAfterGrowth.h>
#include <CreateVoxelSpaceData.h>
///\file CreateHDF5Files.h
/// \brief HDF5 files for simulation results.
///
/// After simulation one HDF5 file is created for simulation results
/// and one HDF5 file for the all trees as XML strings.  An HDF5 file
/// is self-documenting and the groups and datasets can be retrieved
/// following the hierarchical tree-like structure of group and dataset names. 
///
/// HDF5 groups and datasets for simulation results contain the simulation configuration
/// including the program command line, data collected from individual trees each simulation year
/// as well as aggregate stand level data.The 2D and 3D datasets have also attribute
/// names for data columns.
///
/// -# /AllFunctionFiles: Group for functions. Function files are saved as strings.
///    Function files under AllFunctions have suffix *.fun*.
/// -# /AllFunctions: Group for functions. Function files are saved as HDF5 2D tables.
/// -# /AllMetaFiles: Group for Lignum MetaFiles. MetaFiles are saved as strings.
/// -# /AllParameterFiles Group for parameter files. Parameter files are saved *as is* strings.
///   -# /AllParameterFiles/Tree*.txt Tree parameter file(s) with traditional *Tree* prefix .
///   -# /AllParameterFiles/dhlimit.txt Define growth limit for first and second order branches.
/// -# /AllParameters: Group for parameters. Parameter files are saved as HDF5 2D tables.
/// -# /CenterStandData: HDF5 2D table for aggregate center stand data
/// -# /CommandLine: Command line used saved as a string
/// -# /Firmament/Firmament*.txt: Firmament dataset saved as a string. Traditional file prefix *Firmament*. 
/// -# /ForestDeadTreeData: HDF5 3D table for data from dead trees.
/// -# /ForestTreeData: HDF5 3D table for data from growing trees
/// -# /Parameters/TreeParameters: Tree parameters (denoted by relevant Lignum::LGMPD names) as HDF5 1D array
/// -# /StandData: HDF5 2D array for aggregate stand data
/// -# /TreeFunctions: Group for functions in a tree. Function files are save as HDF5 2D tables
/// -# /VoxelSpace/VoxelSpace*.txt: The initial voxel space file. Traditional file prefix *VoxelSpace*.
/// -# /VoxelSpaceSizes/VoxelSpaceSizesData: HDF5 2D data for voxel space expansion during simulation years.
///
/// Trees are saved in an HDF5 file as XML strings and grouped by simulation years. Dataset names are
/// based on unique tree identification tags. Trees are collected by pertinent user defined intervals. For example:
///
/// -# /TreeXML: Root group for XML trees
/// -# /TreeXML/10: Group for trees collected simulation year 10
/// -# /TreeXML/20/Tree_967: XML string for tree simulation year 20 with tree identification tag 967
///
/// The HDF5 files are meant to be analysed with R, Python or some other high level language or tool.
/// R and Python have packages (*rhdf5* and *h5py* respectively) implementing HDF5 API to read, write and examine
/// HDF5 files and file contents.

namespace LignumForest{
  ///\brief After simulation create HDF5 file for simulation results.
  ///
  ///The HDF5 file contains also parameters, functions used as well as the command line
  ///\pre Simulation i.e. the growth loop is done and the results are in `hdf5_data` and in `hdf5_tree_param_data`
  ///\param hdf5fname File name for the HDF5 file
  ///\param hdf5_data Simulation results for each tree and each year, 3D matrix
  ///\param hdf5_tree_param_data Parameter values collected from the trees, 2D matrix
  ///\param argc Command line parameter: number of command line parameters
  ///\param argv Command line parameter: parameters from the command line
  ///\post The HDF5 file `hdf5fname` is closed
  ///\note The string datasets for files used in simulations are collected with wild card search according to customary naming
  /// in the function.
  /// -# *.fun *Function files used*
  /// -# Tree*.txt *Parameters used*
  /// -# MetaFile*.txt *MetaFiles used (MetaFiles has mandatory files needed in simulation)*
  /// -# Firmament*.txt *The Firmament (possibly many) used in simulation*
  /// -# VoxelSpace.txt *The initial Voxel space used*
  /// \deprecated Class CreateHDF5File is in use
  /// \sa CreateHDF5File::createDataSets
  void CreateLignumHDF5File(const string& hdf5fname, const TMatrix3D<double>& hdf5_data, TMatrix2D<double> hdf5_tree_param_data,int argc, char** argv);

  ///\brief Create HDF5 file for simulation data
  ///
  /// After simulation loop collect:
  /// -# Data from each tree, each year
  /// -# Data from dead trees, each year 
  /// -# Stand level aggregate data
  /// -# Center stand aggregate data
  /// -# Meta files used to configure trees
  /// -# Current parameters from trees
  /// -# Current functions from trees
  /// -# All parameters found in simulation directory
  /// -# All functions found in simulation  directory
  /// -# Parameter files as files
  /// -# Function files as files
  /// -# Firmament file as file
  /// -# The initial voxel space file as file
  /// -# VoxelSpace dimensions collected during simulation
  /// -# LignumForest::SEGMENT_LENGTH_LIMIT_FILE
  class CreateHDF5File{
  public:
    /// \brief Create and initialize HDF5 file with its groups.
    /// \param hdf5fname HDF5 file name
    /// \param vsfname Voxel space file name
    /// \param metafile_q MetaFiles used in simulation
    CreateHDF5File(const string& hdf5fname,const string& vsfname, const deque<string>& metafile_q);
    /// Close HDF5 file in destructor
    ~CreateHDF5File();
    /// \brief Create datasets for simulation session configuration:
    /// -# Functions (*.fun files)
    /// -# MetaFile*.txt files
    /// -# Tree*.txt configuration files
    /// -# Firmamanet*.txt files
    /// -# The initial VoxelSpace configuration
    /// -# Max segment length file
    /// -# The command line used
    /// \param argc Number of strings (arguments) in the command line
    /// \param argv Vector of strings in the command line
    /// \attention The quotation characters to protect Glob expressions seems
    /// to be removed by the *copy* algorithm in c++ STL library.
    /// \sa GrowthLoop::parseCommadLine()
    void createConfigurationDataSets(int argc, char** argv);
    /// \brief Create datasets from simulation session results
    /// 
    /// \param gloop Growth loop containing simulation data
    /// Create the following datasests
    /// -# Tree by tree data for each simulation yesr
    /// -# Data from dead trees
    /// -# Stand data
    /// -# Center stand data
    /// -# Tree parameters used
    /// -# Tree functions used
    /// -# Lignum::VoxelSpace dimensions during the simulation 
    template <class T>
    void createDataSets(T& gloop);
    /// Close HDF5 file
    void close();
  private:
    void createGroups();
    LGMHDF5File hdf5_file;///< HDF5 file for Lignum tree and simulations
    string vsfile; ///< Name of the VoxelSpace
    deque<string> metafile_queue; ///< Names for the MetaFiles used
  };

  template <class T>
  void CreateHDF5File::createDataSets(T& gloop)
  {
    //3D dataset for simulation results: each tree, each year (i.e. time is the 3rd dimension)
    const TMatrix3D<double>& hdf5_data = gloop.getHDF5TreeData();
    hdf5_file.createDataSet(LignumForest::TREE_DATA_DATASET_NAME,hdf5_data.rows(),hdf5_data.cols(),hdf5_data.zdim(),hdf5_data);
    hdf5_file.createColumnNames(LignumForest::TREE_DATA_DATASET_NAME,LignumForest::TREE_DATA_COLUMN_ATTRIBUTE_NAME,
				LignumForest::TREE_DATA_COLUMN_NAMES);
    //3D dataset for dead trees: each dead tree appears once the year it has died and removed from the simullation
    const  TMatrix3D<double>& hdf5_dead_tree_data = gloop.getHDF5DeadTreeData();
    hdf5_file.createDataSet(LignumForest::DEAD_TREE_DATA_DATASET_NAME,
			    hdf5_dead_tree_data.rows(),hdf5_dead_tree_data.cols(),hdf5_dead_tree_data.zdim(),hdf5_dead_tree_data);
    hdf5_file.createColumnNames(LignumForest::DEAD_TREE_DATA_DATASET_NAME,LignumForest::DEAD_TREE_DATA_COLUMN_ATTRIBUTE_NAME,
				LignumForest::DEAD_TREE_DATA_COLUMN_NAMES);
    //Stand data
    const TMatrix2D<double>& hdf5_stand_data = gloop.getHDF5StandData();
    hdf5_file.createDataSet(STAND_DATA_DATASET_NAME,hdf5_stand_data.rows(),hdf5_stand_data.cols(),hdf5_stand_data);
    hdf5_file.createColumnNames(STAND_DATA_DATASET_NAME,STAND_DATA_COLUMN_ATTRIBUTE_NAME,STAND_DATA_COLUMN_NAMES);
    //Center stand data
    const TMatrix2D<double>& hdf5_center_stand_data = gloop.getHDF5CenterStandData();
    hdf5_file.createDataSet(CENTER_STAND_DATA_DATASET_NAME,hdf5_center_stand_data.rows(),hdf5_center_stand_data.cols(),
			    hdf5_center_stand_data);
    hdf5_file.createColumnNames(CENTER_STAND_DATA_DATASET_NAME,STAND_DATA_COLUMN_ATTRIBUTE_NAME,STAND_DATA_COLUMN_NAMES);
    //2D dataset for parameters used, each tree
    TMatrix2D<double> hdf5_tree_param_data = gloop.getHDF5TreeParameterData();
    hdf5_file.createDataSet(LignumForest::PGROUP+LignumForest::TREE_PARAMETER_DATASET_NAME,hdf5_tree_param_data.rows(),hdf5_tree_param_data.cols(),
			    hdf5_tree_param_data);
    hdf5_file.createColumnNames(LignumForest::PGROUP+LignumForest::TREE_PARAMETER_DATASET_NAME,LignumForest::TREE_PARAMETER_ATTRIBUTE_NAME,
				LignumForest::TREE_PARAMETER_NAMES);
    //Tree functions
    for (unsigned int i=0; i < FN_V.size();i++){ 
      TMatrix2D<double> hdf5_tree_fn_data = gloop.getHDF5TreeFunctionData(FN_V[i]);
      hdf5_file.createDataSet(TFGROUP+FNA_STR[i],hdf5_tree_fn_data.rows(),hdf5_tree_fn_data.cols(),hdf5_tree_fn_data);
      hdf5_file.createColumnNames(TFGROUP+FNA_STR[i],TREE_FN_ATTRIBUTE_NAME,TREE_FN_COLUMN_NAMES);
    }

    //VoxelSpace dimensions during simulation
    const CreateVoxelSpaceData&  vsdata = gloop.getVoxelSpaceData();
    const TMatrix2D<double>& vsmatrix = vsdata.getData();
    hdf5_file.createDataSet(LignumForest::VOXELSPACESIZESGROUP+LignumForest::VOXELSPACESIZES_DATASET_NAME,vsmatrix.rows(),vsmatrix.cols(),vsmatrix);
    hdf5_file.createColumnNames(LignumForest::VOXELSPACESIZESGROUP+LignumForest::VOXELSPACESIZES_DATASET_NAME,LignumForest::VOXELSPACESIZES_ATTRIBUTE_NAME,
				LignumForest::VS_SIZES_COLUMN_NAMES);
  }
}//namespace LignumForest
#endif
