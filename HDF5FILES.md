# HDF5 files

After simulation one HDF5 file is created for simulation results
and one HDF5 file for trees as XML strings.  An HDF5 file
is a self-documenting tree-like collection of groups and datasets.

### Simulation results
HDF5 file for simulation results contains the simulation configuration
including the program command line, data collected from individual trees each simulation year
as well as aggregate stand level data. In short, all data needed to repeat
the simulation is compiled. The 2D and 3D datasets have also attribute names for data columns 
and are essential for the self-documentation. 

#### Command line
- <I> /CommandLine</I>: Command line string used.

#### Stand data
- <I>/ForestTreeData</I>: HDF5 3D table for data from growing trees
- <I>/ForestDeadTreeData</I>: HDF5 3D table for data from dead trees.
- <I>/StandData</I>: HDF5 2D table for aggregate stand data
- <I>/CenterStandData</I>: HDF5 2D table for center stand data.

#### Functions and parameters
- <I>/AllFunctionFiles</I>: Group for function files with the traditional file suffix <I>.fun</I>. Function files are saved as strings.
- <I>/AllMetaFiles</I>: Group for Lignum MetaFiles. MetaFiles are saved as strings, assuming traditional *MetaFile* file prefix.
- <I>/AllParameterFiles</I>: Group for parameter files. Parameter files are saved as strings.
  - <I>/AllParameterFiles/Tree\*.txt</I>: Tree parameter file(s) assuming the traditional *Tree* file prefix .
  - <I>/AllParameterFiles/dhlimit.txt</I>: The file to define growth limit for the first and second order branches.

#### Firmament
- <I>/Firmament/Firmament\*.txt</I>: Firmament configuration saved as a string assuming the traditional *Firmament* file prefix. 

#### Voxel space
- <I>/VoxelSpace/VoxelSpace\*.txt</I>: The initial voxel space file assuming the traditional *VoxelSpace* file prefix.
- <I>/VoxelSpaceSizes/VoxelSpaceSizesData</I>: HDF5 2D table for voxel space expansion during simulation years.

#### Supplementary data
The supplementary data duplicate their corresponding files saved and can be removed from the HDF5 file later.

- <I>/AllFunctions</I>: Group for functions. Functions are saved as HDF5 2D tables.
- <I>/TreeFunctions</I>: Group for functions in a tree. Function files are save as HDF5 2D tables
- <I>/AllParameters</I>: Group for parameters. Parameter values are saved as HDF5 2D tables.
- <I>/Parameters/TreeParameters</I>: Tree parameters (denoted by relevant Lignum::LGMPD names) as HDF5 1D array.

### Trees
Trees are saved in an HDF5 file as XML strings and grouped by simulation years. The file
name is the file name for simulation results prefixed with *TreesXML_*.

Dataset names for trees are based on unique tree identification tags. Trees are collected 
by user defined intervals in years, for example:

- <I>/TreeXML</I>: Root group for XML trees
- <I>/TreeXML/10</I>: Group for trees collected simulation year 10
- <I>/TreeXML/20/Tree_967</I>: XML string for tree simulation year 20 with tree identification tag 967

See `lignum-forest ` command line for details. 

The HDF5 files are meant to be analysed with R, Python or some other high level language or tool.
R and Python have packages (*rhdf5* and *h5py* respectively) implementing HDF5 API to read, write and examine
HDF5 files and file contents.

## Extract data from HDF5 files
The examples assume simulations are done in the *LignumForest* working directory.

### Simulation results 
To extract simulation configuration data and create the pdf result files
use `ResultAnalysis/ForestPlotAndConfig.R`, for example:
	
	R
	>source('ResultAnalysis/ForestPlotAndConfig.R',chdir=TRUE)
	>ForestPlotAndConfig('SimulationResults.h5,pick=5,GYdata='ResultAnalysis',outdir='ResultsDir')
	>
	
The example creates data files in *ResultsDir* (*outdir=ResultsDir*) using every 5th tree (*pick*=5). 
See the file `ForestPlotAndConfig.R` how it collects and uses several other result analysis files 
as a single batch.

### Trees
Trees are saved as XML strings in user defined intervals. To extract trees use 
`ResultAnalysis/ExtractXML.R`, for example:

	R
	>source('ResultAnalysis/ExtractXML.R')
	>setwd('ResultsDir')
	>ExtractXML('SimulationResults.h5','TreesXML_SimulationResults.h5',c(20,40,60))
	>

The example extracts and creates XML tree files for the shortest, median and the longest trees
for the years 20,40 and 60. See `ExtractXML.R` for details.


