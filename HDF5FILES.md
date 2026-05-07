# HDF5 files
Each simulation run creates three HDF5 files — one for the primary results, a second for tree structures 
stored as XML strings and the third for the voxel space data.

## File names
File names are based on the primary HDF5 result file, with *TreesXML* or *VoxelSpaces* append for tree and
voxel space data.

## Primary simulation results
To ensure reproducibility, all primary simulation data — ranging from the initial configuration
and command line to annual tree and stand-level outputs — is compiled into one HDF5 file. 
Furthermore, the 2D and 3D datasets are self-documenting through the use of descriptive 
attribute names

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
- <I>/VoxelSpaceSizes/VoxelSpaceSizesData</I>: HDF5 2D table for voxel space size expansion during simulation.

Dataset saves voxel space size dynamics. See HDF5 voxel space file for tree, foliage and solar radiation data.

### Supplementary data
The obsolete supplementary data duplicate the content of the files and will be removed from the HDF5 file structure.

- <I>/AllFunctions</I>: Group for functions. Functions are saved as HDF5 2D tables.
- <I>/TreeFunctions</I>: Group for functions in a tree. Function files are save as HDF5 2D tables.
- <I>/AllParameters</I>: Group for parameters. Parameter values are saved as HDF5 2D tables.
- <I>/Parameters/TreeParameters</I>: Tree parameters (denoted by relevant Lignum::LGMPD names) as HDF5 1D array.

## Tree structures
Trees are saved in a separate HDF5 file as XML strings and grouped by simulation years. 
The designated file name is the file name for simulation results prefixed with *TreesXML_*.

Dataset names for trees are based on unique tree identification tags. Trees are collected 
by user defined intervals for certain years, for example:

- <I>/TreeXML</I>: Main group for the trees.
- <I>/TreeXML/10</I>: Group for trees collected for the simulation year 10.
- <I>/TreeXML/20/Tree_967</I>: XML string for the tree collected for the simulation year 20 with the tree identification tag 967.

## Voxel space data
*VoxelSpace* data sets are collected at regular intervals as 4D matrices, where the fourth dimension represents 
the tree, foliage and solar radiation data recorded from each voxel. 

## Retrieve data from HDF5 files

### Simulation analyses
The `R` scripts in `ResultAnalysis/ForestPlotAndConfig.R` analyze LignumForest simulation. 
*SimulationResults.h5* serves as an example containing the primary simulation data:
	
	R #Start R
	>source('ResultAnalysis/ForestPlotAndConfig.R',chdir=TRUE)
	>ForestPlotAndConfig('SimulationResults.h5',pick=5,GYdata='ResultAnalysis/',outdir='ResultsDir')
	>
	
`ForestPlotAndConfig` generates PDF plots and retrieves simulation configurations. These files, along with 
the simulation’s HDF5 data, are moved to the *outdir* directory. To optimize processing, the *pick* parameter 
evaluates every fifth tree, while *GYdata* specifies the location of the growth and yield tables required 
for the analysis.

### Trees
Continuing the example Trees are saved as XML strings in the file *TreesXML_SimulationResults.h5*. 
Exctract them with `ResultAnalysis/ExtractXML.R` script:
 
	>source('ResultAnalysis/ExtractXML.R')
	>setwd('ResultsDir') #Move to 'outdir' directory set in ForestPlotAndConfig
	>ExtractXML('SimulationResults.h5','TreesXML_SimulationResults.h5',c(20,40,60))
	>

`ExtractXML` generates XML files for the shortest, median, and longest trees for the years 20, 40, and 60.
 For larger specimen sets, use LignumVTK.

### Visualization
Use the LignumVTK project to visualize trees and voxel space data with ParaView.
