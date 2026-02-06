# HDF5 files
After simulation one HDF5 file is created for simulation results
and one HDF5 file for trees as XML strings.  An HDF5 file
is a self-documenting tree-like hierarchical collection of datasets 
and groups of datasets.

The HDF5 files are meant to be analysed with R, Python or some other 
high level language or tool. R and Python have packages *rhdf5* and *h5py* 
respectively implementing HDF5 API to read, write and examine HDF5 files.

## Simulation data
The HDF5 file for the simulation data contains the simulation configuration
including the program command line, data collected from individual trees 
each simulation year as well as aggregate stand level data. In short, simulation 
data collected and all data needed to repeat the simulation is compiled into a single HDF5 file. 
The 2D and 3D datasets have also attribute names for data columns for the self-documentation. 

### Command line
- <I>/CommandLine</I>: Command line string used.

### Stand data
- <I>/ForestTreeData</I>: HDF5 3D table for data from growing trees.
- <I>/ForestDeadTreeData</I>: HDF5 3D table for data from dead trees. 
  A dead tree appears once, in the year when it is removed from the simulation.
- <I>/StandData</I>: HDF5 2D table for aggregate stand data.
- <I>/CenterStandData</I>: HDF5 2D table for center stand data.

### Functions and parameters
- <I>/AllFunctionFiles</I>: Group for function files with the traditional file suffix <I>.fun</I>. Function files are saved as strings.
- <I>/AllMetaFiles</I>: Group for Lignum MetaFiles. MetaFiles are saved as strings assuming the traditional *MetaFile* file prefix.
- <I>/AllParameterFiles</I>: Group for parameter files. Parameter files are saved as strings.
  - <I>/AllParameterFiles/Tree\*.txt</I>: Tree parameter files assuming the traditional *Tree* file prefix .
  - <I>/AllParameterFiles/dhlimit.txt</I>: The file to define growth limit for the first and second order branches.

### Firmament
- <I>/Firmament/Firmament\*.txt</I>: Firmament configuration saved as a string assuming the traditional *Firmament* file prefix. 

### Voxel space
- <I>/VoxelSpace/VoxelSpace\*.txt</I>: The initial voxel space assuming the traditional *VoxelSpace* file prefix.
- <I>/VoxelSpaceSizes/VoxelSpaceSizesData</I>: HDF5 2D table for voxel space expansion during simulation.

### Supplementary data
The obsolete supplementary data duplicate the content of their corresponding files saved 
and can be removed from the HDF5 file structure some time in the future.

- <I>/AllFunctions</I>: Group for functions. Functions are saved as HDF5 2D tables.
- <I>/TreeFunctions</I>: Group for functions in a tree. Function files are save as HDF5 2D tables.
- <I>/AllParameters</I>: Group for parameters. Parameter values are saved as HDF5 2D tables.
- <I>/Parameters/TreeParameters</I>: Tree parameters (denoted by relevant Lignum::LGMPD names) as HDF5 1D array.

### Trees
Trees are saved in a separate HDF5 file as XML strings and grouped by simulation years. 
The designated file name is the file name for simulation results prefixed with *TreesXML_*.

Dataset names for trees are based on unique tree identification tags. Trees are collected 
by user defined intervals for certain years, for example:

- <I>/TreeXML</I>: Main group for the trees.
- <I>/TreeXML/10</I>: Group for trees collected for the simulation year 10.
- <I>/TreeXML/20/Tree_967</I>: XML string for the tree collected for the simulation year 20 with the tree identification tag 967.

## Retrieve data from HDF5 files
The example for simulation analyses and trees assume *LignumForest* working directory.

### Simulation analyses
`ResultAnalysis/ForestPlotAndConfig.R` is a collection of `R` scripts 
for LignumForest simulation analyses. Suppose the HDF5 file *SimulationResults.h5* 
contains the data from a simulation:
	
	R #Start R
	>source('ResultAnalysis/ForestPlotAndConfig.R',chdir=TRUE)
	>ForestPlotAndConfig('SimulationResults.h5',pick=5,GYdata='ResultAnalysis/',outdir='ResultsDir')
	>
	
`ForestPlotAndConfig` creates pdf files for analyses results and retrieves simulation configuration. 
These as well as the HDF5 files from the simulation are moved to *ResultsDir*. The *pick* parameter 
selects every 5th tree for evaluation. The *GYdata* parameter points to the directory with 
predefined growth and yield tables used in analyses. 

### Trees
Trees are saved as XML strings. Continuing the example the trees are in the file
*TreesXML_SimulationResults.h5*. To extract the trees use `ResultAnalysis/ExtractXML.R`:
 
	>source('ResultAnalysis/ExtractXML.R')
	>setwd('ResultsDir') #Move to 'outdir' directory set in ForestPlotAndConfig
	>ExtractXML('SimulationResults.h5','TreesXML_SimulationResults.h5',c(20,40,60))
	>

`ExtractXML` creates XML files for the shortest, the median and the longest tree
for the years 20,40 and 60. If more than three specimen trees for each selected 
year are needed for inspection see the LignumVTK project to select trees for visualization.

