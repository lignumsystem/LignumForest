#include <CreateHDF5Files.h>
///\file CreateHDF5Files.cc
namespace LignumForest{
  void CreateLignumHDF5File(const string& hdf5fname, const TMatrix3D<double>& hdf5_data, TMatrix2D<double> hdf5_tree_param_data,int argc, char** argv)
  {
    //Create and open HDF5 file for write
    LGMHDF5File hdf5_file(hdf5fname);
    //Create dataset groups
    //Parameters
    hdf5_file.createGroup(LignumForest::PGROUP);
    //All functions (*.fun) 
    hdf5_file.createGroup(LignumForest::AFGROUP);
    hdf5_file.createGroup(LignumForest::PFILEGROUP);
    hdf5_file.createGroup(LignumForest::ALLMETAFILEGROUP);
    hdf5_file.createGroup(LignumForest::ALLPARAMFILEGROUP);
    hdf5_file.createGroup(LignumForest::ALLFNFILEGROUP);
    hdf5_file.createGroup(LignumForest::FIRMAMENTGROUP);
    hdf5_file.createGroup(LignumForest::VOXELSPACEGROUP);
    //Data sets
    //3D dataset for simulation results: each tree, each year (i.e. time is the 3rd dimension)
    hdf5_file.createDataSet(LignumForest::TREE_DATA_DATASET_NAME,hdf5_data.rows(),hdf5_data.cols(),hdf5_data.zdim(),hdf5_data);
    hdf5_file.createColumnNames(LignumForest::TREE_DATA_DATASET_NAME,LignumForest::TREE_DATA_COLUMN_ATTRIBUTE_NAME,LignumForest::TREE_DATA_COLUMN_NAMES);
    //2D dataset for parameters used, each tree
    hdf5_file.createDataSet(LignumForest::PGROUP+LignumForest::TREE_PARAMETER_DATASET_NAME,hdf5_tree_param_data.rows(),hdf5_tree_param_data.cols(),
			    hdf5_tree_param_data);
    hdf5_file.createColumnNames(LignumForest::PGROUP+LignumForest::TREE_PARAMETER_DATASET_NAME,LignumForest::TREE_PARAMETER_ATTRIBUTE_NAME,
				LignumForest::TREE_PARAMETER_NAMES);
    //All functions from files in a simulation directory as 2D dataset
    hdf5_file.createFnDataSetsFromDir("*.fun",LignumForest::AFGROUP,LignumForest::TREE_FN_ATTRIBUTE_NAME,LignumForest::TREE_FN_COLUMN_NAMES);
    //All Tree parameters in a simulation directory  as 2D dataset
    hdf5_file.createParameterDataSetsFromDir("Tree*.txt",LignumForest::PFILEGROUP,LignumForest::TREE_PARAMETER_FILE_ATTRIBUTE_NAME,
					     LignumForest::TREE_PARAMETER_COLUMN_NAMES);
    //String datasets for files used
    //Metafiles
    hdf5_file.createFileDataSetsFromDir("MetaFile*.txt",LignumForest::ALLMETAFILEGROUP);
    //Tree parameters
    hdf5_file.createFileDataSetsFromDir("Tree*.txt",LignumForest::ALLPARAMFILEGROUP);
    //All functions
    hdf5_file.createFileDataSetsFromDir("*.fun",LignumForest::ALLFNFILEGROUP);
    //The Firmament used
    hdf5_file.createFileDataSetsFromDir("Firmament*.txt",LignumForest::FIRMAMENTGROUP);
    //The Voxel space used
    hdf5_file.createFileDataSetsFromDir("VoxelSpace.txt",LignumForest::VOXELSPACEGROUP);
    //Command line
    vector<string> c_vec;
    std::copy( argv, argv+argc,back_inserter(c_vec));
    ostringstream cline;
    copy(c_vec.begin(),c_vec.end(),ostream_iterator<string>(cline, " "));
    hdf5_file.createDataSet(LignumForest::COMMAND_LINE_DATASET_NAME,cline.str());
    //Close HDF5 file
    hdf5_file.close();
  }
  
}	     
