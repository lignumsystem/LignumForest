# LignumForest
LignumForest simulates a growing Scots pine tree community using individual LIGNUM trees. 
There are two shell scripts to run the `lignum-forest` program:
	
    run-lignum-forest.sh
	run-lignum-forest.slurm #SLURM workload manager additions

Run `lignum-forest` without parameters to see all command line options."

## CMake Makefile build system
See the [README](https://github.com/lignumsystem/lignum-core/blob/master/README.md) for
lignum-core. The `make install`command will copy `lignum-forest` to LignumForest directory.

## CMake Xcode build system

### Create Xcode project 
Create the Xcode project file:
	
	cd LignumForest
    mkdir xcode
    cd xcode
    cmake .. -G Xcode

Open the project file in Xcode:
     
	 open LignumForest.xcodeproj

### Build lignum-forest binary
First set Scheme to `lignum-forest`:

	Xcode -> Product (in the menu bar) -> Scheme 

Build the `lignum-forest` binary in  Xcode for Testing. It will appear
in the *xcode/Debug*  directory:

	Xcode -> Product (in the menu bar) -> Build For -> Running/Testing/Profiling

Xcode tracks source file dependencies during the build process. To install 
`lignum-forest` to LignumForest working directory select the `install` Scheme
and build:

	Xcode -> Product (in the menu bar) -> Build 
	
### Debugging the program
Copy necessary function and parameter files (*.fun* and *.txt* suffixes)
to *xcode/Debug*  to the `lignum-forest` binary location. `lignum-forest` assumes 
that the configuration files are found in the working directory. 

Set Scheme  to `lignum-forest`:

	Xcode -> Product (in the menu bar) -> Scheme 

Set command  line parameters for  `lignum-forest` in Xcode:

	Xcode -> Product (in the menu  bar) -> Scheme ->  Edit Scheme -> Arguments.

Divide the lengthy command line into practical parts for debugging from `Arguments -> '+'`.

Set the breakpoints in source files. To debug the program:

	Xcode -> Product (in the menu bar) -> Run

Alternatively, load the `lignum-forest` binary to Xcode from the LignumForest directory:

	Xcode -> Debug (in the menu bar) -> Debug Executable

## Software documentation
Run `doxygen` in the LignumForest directory:
    
    doxygen Doxyfile 2> errors.txt
     
To read html version of the document type:

    open DoxygenDoc/html/index.html
    
Go to latex subdirectory and use make:

    cd DoxygenDoc/latex
    make all
	open refman.pdf
    
The final LaTeX documentation will be provided as *refman.pdf*. 

The main growth loop for LignumForest is implemented in *lignum-forest.cc*.
Simulation results are saved in [HDF5 files](HDF5FILES.md). The introductionary presentation
will appear in [GENERAL_DESCRIPTION](GENERAL_DESCRIPTION.md).

> [!TIP]
> Doxyfile default values for graphs may result too large uninformative figures.
> Reduce values, for example DOT_GRAPH_MAX_NODES = 10 and MAX_DOT_GRAPH_DEPTH = 3,
> for more concise graphs and network diagrams easier to understand figures.

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
