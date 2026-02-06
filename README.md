# LignumForest
LignumForest is a project for simulating a growing tree community using individual LIGNUM trees. 
The main growth loop for LignumForest is implemented in *lignum-forest.cc*.
Simulation results are saved in [HDF5 files](HDF5FILES.md).

## CMake build system

### CMake for Makefile build system
Create Makefile build system for the *lignum-core* system. 
See the [README](https://github.com/lignumsystem/lignum-core/blob/master/README.md)
and the *CMakeLists.txt* file in *lignum-core* for details.

To create Makefile build system for LignumForest to compile 
`lignum-forest` binary for debugging:

	git clone https://github.com/lignumsystem/LignumForest.git #Download the software
    cd LignumForest
    mkdir debug 
    cd  debug
    cmake .. -DCMAKE_BUILD_TYPE=Debug
    make install 

Makefile build system for release (optimised code, no debug information):

    cd LignumForest
    mkdir release
    cd release
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make install

The `make install`command will copy `lignum-forest` to LignumForest directory
where there are two shell scripts to run the program:
	
    run-lignum-forest.sh
	run-lignum-forest.slurm
	
Command line options and their  short documentation can be obtained by
running `./lignum-forest`  without any  command line parameters.  See also
LignumForest::GrowthLoop::usage().

> [!IMPORTANT]
> It is important to type `make install` to move `lignum-forest` to
> directory above to be used by the scripts to run simulatations.
> Typing just `make` the `lignum-forest` program remains in the compilation directory.

To recompile `lignum-forest` type:

	make clean
	make install
	
> [!IMPORTANT]
> CMake tracks by default file changes only in the current project (LignumForest in this case). 
> To let CMake  follow all file dependencies correctly `make clean` is mandatory before recompilation. 
> After `make clean` CMake will have correct build tree from previous software  build.

> [!NOTE]
> To remove all CMake  configurations and compilation work just
> remove the build  directory (i.e. *debug*,  *release* or *xcode*)
> and recreate it.

### CMake for Xcode

#### Create Xcode project 
Create the Xcode project file:
	
	cd LignumForest
    mkdir xcode
    cd xcode
    cmake .. -G Xcode

Open the project file in Xcode:
     
	 open LignumForest.xcodeproj

#### Build lignum-forest binary
First set Scheme to `lignum-forest`:

	Xcode -> Product (in the menu bar) -> Scheme 

Build the `lignum-forest` binary in  Xcode for Testing. It will appear
in the *xcode/Debug*  directory:

	Xcode -> Product (in the menu bar) -> Build For -> Running/Testing/Profiling

Xcode tracks source file dependencies during the build process. To install 
`lignum-forest` to LignumForest working directory select the `install` Scheme
and build:

	Xcode -> Product (in the menu bar) -> Build 
	
#### Debugging the program
Copy necessary function files and parameter files (*.fun* and *.txt* suffixes)
to *xcode/Debug*  where the `lignum-forest` binary is  located. 
`lignum-forest` assumes that the configuration files 
are found in the working directory. 

Make sure that Scheme is set to `lignum-forest`:

	Xcode -> Product (in the menu bar) -> Scheme 

Set command  line parameters for  `lignum-forest` in Xcode:

	Xcode -> Product (in the menu  bar) -> Scheme ->  Edit Scheme -> Arguments.

Divide the lengthy command line into practical parts for debugging from `Arguments -> '+'`.

Set the breakpoints in source files. To debug the program:

	Xcode -> Product (in the menu bar) -> Run

Alternatively, load the `lignum-forest` binary to Xcode from the LignumForest directory:

	Xcode -> Debug (in the menu bar) -> Debug Executable

### LignumForest project dependency graph
CMake can generate Graphviz file representing library and binary dependencies in the project 
using the Dot language grammar. After that create an image file using `dot`. 
In the *LignumForest* directory create the *graphviz* build directory:
	
	mkdir graphviz
	cd graphviz
	cmake ..   --graphviz=LignumForest.dot
	dot -Tpdf -Kneato -Goverlap=prism  LignumForest.dot  -o  LignumForest.pdf
	
The file *LignumForest.pdf* contains the visual presentation of the Graphviz file. 
The option `-T` supports many well known image file formats.

## Software documentation
The introductionary presentation will appear in [GENERAL_DESCRIPTION](GENERAL_DESCRIPTION.md).

The Reference Guide for the LignumForest will be based on Doxygen typesetting and other information
available for the software. To generate the documentation run `doxygen` in the LignumForest directory:
    
    doxygen Doxyfile 2> errors.txt
     
Doxyfile is the configuration file for `doxygen`. The documentation will appear in DoxygenDoc directory. 
Errors and warnings will appear in *errors.txt*. To read html version of the document type (macOS):

    open DoxygenDoc/html/index.html
    
To generate LaTeX version go to latex subdirectory and use make:

    cd DoxygenDoc/latex
    make all
	open refman.pdf
    
The result will be *refman.pdf* that can be opened with a pdf reader.

> [!TIP]
> Doxyfile default values for graphs may result too large uninformative figures.
> Reduce for example DOT_GRAPH_MAX_NODES = 10 and MAX_DOT_GRAPH_DEPTH = 3 
> for more concise graphs and network diagrams easier to understand figures.

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
