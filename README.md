# LignumForest
LignumForest simulates a growing Scots pine tree community using individual LIGNUM trees. 
There are two shell scripts to run the `lignum-forest` program:
	
    run-lignum-forest.sh
	run-lignum-forest.slurm #Slurm workload manager additions

The Slurm[^slurm] version includes native workload management commands.

`lignum-forest` looks for configuration files in the working directory. 
Run `lignum-forest` without parameters to see all command line options."

The *lignum-forest.cc* file handles the main growth loop. Simulation data is saved in HDF5 format,
and an introductory overview will be in GENERAL_DESCRIPTION.

## Prerequisites
Download *lignum-core* and *LignumForest* from GitHub. Install the dependencies specified in the 
lignum-core [README](https://github.com/lignumsystem/lignum-core/blob/master/README.md). 

## CMake Makefile build system
Running `make install` compiles `lignum-forest` and copies it into the LignumForest working directory.

## CMake Xcode build system
### Create Xcode project 
Create the Xcode project file and open the project file in Xcode:
     
	 open LignumForest.xcodeproj

The CMakeLists.txt file defines the Xcode project's targets, sources, and build settings.

### Build lignum-forest
Set Scheme to `lignum-forest`:

	Xcode -> Product -> Scheme -> lignum-forest

Build the `lignum-forest` binary in  Xcode for Testing. It will appear
in the *xcode/Debug*  directory:

	Xcode -> Product -> Build For -> Testing
	
### Debugging the program
`lignum-forest` looks for configuration files in the working directory. Copy the required
function and parameter files (*.fun and *.txt) to the *xcode/Debug* folder where the 
`lignum-forest` binary is located.

Set Scheme  to `lignum-forest`:

	Xcode -> Product -> Scheme -> lignum-forest
	
Break down the long command line for easier debugging:

	Xcode -> Product -> Scheme ->  Edit Scheme -> Arguments -> '+'.

Set the breakpoints in source files. To debug the program:

	Xcode -> Product -> Run

#### Install lignum-forest
Set Scheme to `install` and build:

	Xcode -> Product -> Build 
	
## Software documentation
Run `doxygen` in the LignumForest directory to generate the html document:
    
    doxygen Doxyfile 2> errors.txt
    open DoxygenDoc/html/index.html
	
Run the Doxygen-generated Makefile to compile the LaTeX output into *refman.pdf*:

    cd DoxygenDoc/latex
    make all
	open refman.pdf

> [!TIP]
> The Doxyfile defaults can produce overly large, uninformative figures. To create more concise
> graphs and network diagrams, try reducing the limit values — for example, set DOT_GRAPH_MAX_NODES = 10
> and MAX_DOT_GRAPH_DEPTH = 3.

### LignumForest project dependency graph
CMake can generate Graphviz files representing project dependencies using the Dot language. 
You can then convert these into images using the `dot` command. To start, create a 
graphviz build directory inside the LignumForest folder:
	
	mkdir graphviz
	cd graphviz
	cmake ..   --graphviz=LignumForest.dot
	dot -Tpdf -Kneato -Goverlap=prism  LignumForest.dot  -o  LignumForest.pdf
	
The option `-T` supports many well known image file formats. The final figure is
in LignumForest.pdf.

## Litterature to cite
The LIGNUM conifer trees and other components of this project have been used for the calculations 
in the following publications:

- R. Sievänen, J. Perttunen, E. Nikinmaa, and P. Kaitaniemi. Toward extension of a single tree functional-structural model 
of Scots pine to stand level: effect of the canopy of randomly distributed, identical trees on development of tree structure. 
Functional Plant Biology, 35(9/10):964–975, 2008.
- R. Sievänen, P. Raumonen, J. Perttunen, E. Nikinmaa, and P. Kaitaniemi. A study of crown development mechanisms 
using a shoot-based tree model and segmented terrestrial laser scanning data. 
Annals of Botany, 122(3):423–434, 2018.

To refer to the core model of LIGNUM in general use:

- Perttunen J, Sievänen R, Nikinmaa R, Salminen H, Saarenmaa H, Väkevä J. 1996. LIGNUM: a tree model based on simple structural units. Annals of Botany 77: 87–98.
- Perttunen J, Sievänen R, Nikinmaa E. 1998. LIGNUM: a model combining the structure and the functioning of trees. Ecological Modelling 108: 189–198.

[^slurm]: Simple Linux Utility for Resource Management
