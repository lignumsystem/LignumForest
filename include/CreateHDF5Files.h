#ifndef CREATEHDF5FILES_H
#define CREATEHDF5FILES_H
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
  void CreateLignumHDF5File(const string& hdf5fname, const TMatrix3D<double>& hdf5_data, TMatrix2D<double> hdf5_tree_param_data,int argc, char** argv);
}

#endif
