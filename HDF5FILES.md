# HDF5 files
## HDF5 files for simulation results.

After simulation one HDF5 file is created for simulation results
and one HDF5 file for the all trees as XML strings.  An HDF5 file
is self-documenting and the groups and datasets can be retrieved
following the hierarchical tree-like structure of group and dataset names. 

HDF5 groups and datasets for simulation results contain the simulation configuration
including the program command line, data collected from individual trees each simulation year
as well as aggregate stand level data.The 2D and 3D datasets have also attribute
names for data columns.

- /AllFunctionFiles: *Group for function files*. Function files are saved as strings.
   Function files under AllFunctions have suffix *.fun*.
- /AllFunctions: *Group for functions*. Functions are saved as HDF5 2D tables.
- /AllMetaFiles: *Group for Lignum MetaFiles*. MetaFiles are saved as strings, traditionally *MetaFile* prefix.
- /AllParameterFiles *Group for parameter files*. Parameter files are saved *as is* strings.
  - /AllParameterFiles/Tree\*.txt Tree parameter file(s) with traditional *Tree* prefix .
  - /AllParameterFiles/dhlimit.txt Define growth limit for first and second order branches.
- /AllParameters: *Group for parameters*. Parameter files are saved as HDF5 2D tables.
- /CenterStandData: *HDF5 2D table* for aggregate center stand data
- /CommandLine: *Command line* used saved as a string
- /Firmament/Firmament\*.txt: *Firmament configuration* saved as a string. Traditional file prefix *Firmament*. 
- /ForestDeadTreeData: *HDF5 3D table* for data from dead trees.
- /ForestTreeData: *HDF5 3D table* for data from growing trees
- /Parameters/TreeParameters: *Tree parameters* (denoted by relevant Lignum::LGMPD names) as HDF5 1D array
- /StandData: *HDF5 2D array* for aggregate stand data
- /TreeFunctions: *Group for functions* in a tree. Function files are save as HDF5 2D tables
- /VoxelSpace/VoxelSpace\*.txt: *The initial voxel space file*. Traditional file prefix *VoxelSpace*.
- /VoxelSpaceSizes/VoxelSpaceSizesData: *HDF5 2D data for voxel space expansion* during simulation years.

Trees are saved in a different HDF5 file as XML strings and grouped by simulation years. The file
name has the prefix *TreesXML_* appended with the file name for the simulation results. See `lignum-forest `
command line for details. Dataset names are based on unique tree identification tags. 
Trees are collected by user defined intervals in years. For example:

- /TreeXML: Root group for XML trees
- /TreeXML/10: Group for trees collected simulation year 10
- /TreeXML/20/Tree_967: XML string for tree simulation year 20 with tree identification tag 967

The HDF5 files are meant to be analysed with R, Python or some other high level language or tool.
R and Python have packages (*rhdf5* and *h5py* respectively) implementing HDF5 API to read, write and examine
HDF5 files and file contents.

## Extract data from HDF5 files

### Simulation results 
To extract simulation configuration data and create the pdf result files
use `ResultAnalysis/ForestPlotAndConfig.R`, for example:
	
	R
	>source('path/to/ResultAnalysis/ForestPlotAndConfig.R',chdir=TRUE')
	>ForestPlotAndConfig('path/to/SimulationResults.h5,pick=5,GYdata='path/to/ResultAnalysis',outdir='ResultsDir')
	>
	
The example creates data files in *ResultsDir* (*outdir=ResultsDir*) using every 5th tree (*pick*=5). 
See the file `ForestPlotAndConfig.R` how it collects and use several other result analysis files 
as a single batch.

### Trees
Trees are saved as XML strings in user defined intervals. To extract trees use 
`ResultAnalysis/ExtractXML.R`, for example:

	R
	>source('ResultAnalysis/ExtractXML.R')
	>ExtractXML('SimulationResults.h5','TreesXML_SimulationResults.h5',c(20,40,60))
	>

The example extracts and creates XML tree files for the years 20,40 and 60. 
See `ExtractXML.R` for details.


