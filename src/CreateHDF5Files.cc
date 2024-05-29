#include <CreateHDF5Files.h>
///\file CreateHDF5Files.cc
namespace LignumForest{
  CreateHDF5File::CreateHDF5File(const string& hdf5fname,const string& vsfilename)
    :hdf5_file(hdf5fname),vsfile(vsfilename)
  {;
    //Create dataset groups
    createGroups();
  }

  CreateHDF5File::~CreateHDF5File()
  {
    hdf5_file.close();
  }
  
  void CreateHDF5File::createGroups()
  {
    //Create dataset groups
    //Parameters
    hdf5_file.createGroup(LignumForest::PGROUP);
    //Functions in a tree
    hdf5_file.createGroup(TFGROUP);
    //All files for functions and parameters.
    hdf5_file.createGroup(LignumForest::AFGROUP);
    hdf5_file.createGroup(LignumForest::PFILEGROUP);
    //File datasets for MetaFiles, parameters, functions, firmament, voxel space 
    hdf5_file.createGroup(LignumForest::ALLMETAFILEGROUP);
    hdf5_file.createGroup(LignumForest::ALLPARAMFILEGROUP);
    hdf5_file.createGroup(LignumForest::ALLFNFILEGROUP);
    hdf5_file.createGroup(LignumForest::FIRMAMENTGROUP);
    hdf5_file.createGroup(LignumForest::VOXELSPACEGROUP);
  }

  void CreateHDF5File::createConfigurationDataSets(int argc, char** argv)
  {
    //String MetaFile datasets for configuration files used
    //Metafiles
    hdf5_file.createFileDataSetsFromDir("{MetaFile,MetaFile[0-9]}.txt",LignumForest::ALLMETAFILEGROUP);
    //Tree parameters
    //All functions from files in a simulation directory as 2D dataset
    hdf5_file.createFnDataSetsFromDir("*.fun",LignumForest::AFGROUP,LignumForest::TREE_FN_ATTRIBUTE_NAME,LignumForest::TREE_FN_COLUMN_NAMES);
    //All Tree parameters in a simulation directory  as 2D dataset
    hdf5_file.createParameterDataSetsFromDir("{Tree,Tree[0-9]}.txt",LignumForest::PFILEGROUP,LignumForest::TREE_PARAMETER_FILE_ATTRIBUTE_NAME,
					     LignumForest::TREE_PARAMETER_COLUMN_NAMES);
    ///String datasets for parameters, functions, Firmament, initial voxel space
    hdf5_file.createFileDataSetsFromDir("{Tree,Tree[0-9]}.txt",LignumForest::ALLPARAMFILEGROUP);
    //All functions
    hdf5_file.createFileDataSetsFromDir("*.fun",LignumForest::ALLFNFILEGROUP);
    //The Firmament used
    hdf5_file.createFileDataSetsFromDir("{Firmament,Firmament[0-9]}.txt",LignumForest::FIRMAMENTGROUP);
    //The initial Voxel space used
    hdf5_file.createFileDataSetsFromDir(vsfile,LignumForest::VOXELSPACEGROUP);
    //Command line
    vector<string> c_vec;
    std::copy( argv, argv+argc,back_inserter(c_vec));
    ostringstream cline;
    copy(c_vec.begin(),c_vec.end(),ostream_iterator<string>(cline, " "));
    hdf5_file.createDataSet(LignumForest::COMMAND_LINE_DATASET_NAME,cline.str());
  }
  
  void CreateHDF5File::close()
  {
    hdf5_file.close();
  }
  
  
  ///Function to create HDF5 groups and datasets (not part of the class CreateHDF5File).
  void CreateLignumHDF5File(const string& hdf5fname, const TMatrix3D<double>& hdf5_data, TMatrix2D<double> hdf5_tree_param_data, int argc, char** argv)
  {
    //Create and open HDF5 file for write
    LGMHDF5File hdf5_file(hdf5fname);
    //Create dataset groups
    //Parameters
    hdf5_file.createGroup(LignumForest::PGROUP);
    //Functions in a tree
    hdf5_file.createGroup(TFGROUP);
    //All files for functions and prameters.
    hdf5_file.createGroup(LignumForest::AFGROUP);
    hdf5_file.createGroup(LignumForest::PFILEGROUP);
    //File datasets for MetaFiles, parameters, functions, firmament, voxel space 
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
    hdf5_file.createParameterDataSetsFromDir("{Tree,Tree[0-9]}.txt",LignumForest::PFILEGROUP,LignumForest::TREE_PARAMETER_FILE_ATTRIBUTE_NAME,
					     LignumForest::TREE_PARAMETER_COLUMN_NAMES);
    //String datasets for files used
    //Metafiles
    hdf5_file.createFileDataSetsFromDir("{MetaFile,MetaFile[0-9]}.txt",LignumForest::ALLMETAFILEGROUP);
    //Tree parameters
    hdf5_file.createFileDataSetsFromDir("{Tree,Tree[0-9]}.txt",LignumForest::ALLPARAMFILEGROUP);
    //All functions
    hdf5_file.createFileDataSetsFromDir("*.fun",LignumForest::ALLFNFILEGROUP);
    //The Firmament used
    hdf5_file.createFileDataSetsFromDir("{Firmament,Firmament[0-9]}.txt",LignumForest::FIRMAMENTGROUP);
    //The Voxel space used
    hdf5_file.createFileDataSetsFromDir("{VoxelSpace,VoxelSpace[0-9]}.txt",LignumForest::VOXELSPACEGROUP);
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
