#ifndef CREATEHDF5FILES_H
#define CREATEHDF5FILES_H
#include <deque>
#include <LGMHDF5File.h>
#include <TreeDataAfterGrowth.h>

///\file CreateHDF5Files.h
namespace LignumForest{
  ///After simulation create HDF5 file for simulation results.
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
  /// -# *.fun Function files used
  /// -# Tree*.txt Parameters used
  /// -# MetaFile*.txt MetaFiles used (MetaFiles tells mandatory files needed in simulation)
  /// -# Firmament*.txt The Firmament (possibly many) used in simulation
  /// -# VoxelSpace.txt The initial Voxel space used
  /// \todo It might be the best to collect, and easiest to see and retrieve afterwards, string datasets only for files used (functions,
  ///       parameters, Metafiles, Firmament etc.
  /// \todo Save snapshots of the contents of the Voxel space as 3D matrix (matrices) for further study
  /// \note For LignumForest see the class LignumForest::CreateHDF5File.
  void CreateLignumHDF5File(const string& hdf5fname, const TMatrix3D<double>& hdf5_data, TMatrix2D<double> hdf5_tree_param_data,int argc, char** argv);

  ///\brief Create HDF5 file for simulation data
  ///
  /// After simulation loop collect:
  /// -# Data from each tree, each year
  /// -# Stand level aggregate data
  /// -# Center stand aggrgate data
  /// -# Meta files used to configure trees
  /// -# Current parameters from trees
  /// -# Current functions from trees
  /// -# All parameters found in simulation directory
  /// -# All functions found in simulation  directory
  /// -# Parameter files as files
  /// -# Function files as files
  /// -# Firmament file as file
  /// -# The initial voxel space file as file
  class CreateHDF5File{
  public:
    /// Create and initialize HDF5 file with its groups.
    /// \param hdf5fname HDF5 file name
    /// \param vsfname Voxel space file name
    /// \param metafile_q MetaFiles used in simulation
    CreateHDF5File(const string& hdf5fname,const string& vsfname, const deque<string>& metafile_q);
    /// Close HDF5 file in destructor
    ~CreateHDF5File();
    /// Create datasets for simulation session configuration
    void createConfigurationDataSets(int argc, char** argv);
    /// Create datasets from simulation session results
    /// \param gloop Growth loop containing simulation data
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
    TMatrix3D<double>& hdf5_data = gloop.getHDF5TreeData();
    hdf5_file.createDataSet(LignumForest::TREE_DATA_DATASET_NAME,hdf5_data.rows(),hdf5_data.cols(),hdf5_data.zdim(),hdf5_data);
    hdf5_file.createColumnNames(LignumForest::TREE_DATA_DATASET_NAME,LignumForest::TREE_DATA_COLUMN_ATTRIBUTE_NAME,
				LignumForest::TREE_DATA_COLUMN_NAMES);
    //Stand data
    TMatrix2D<double>& hdf5_stand_data = gloop.getHDF5StandData();
    hdf5_file.createDataSet(STAND_DATA_DATASET_NAME,hdf5_stand_data.rows(),hdf5_stand_data.cols(),hdf5_stand_data);
    hdf5_file.createColumnNames(STAND_DATA_DATASET_NAME,STAND_DATA_COLUMN_ATTRIBUTE_NAME,STAND_DATA_COLUMN_NAMES);
    //Center stand data
    TMatrix2D<double>& hdf5_center_stand_data = gloop.getHDF5CenterStandData();
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
  }
}//namespace LignumForest
#endif
