## LignumForest
LignumForest is a project for simulating a growing tree community using individual LIGNUM trees. The LIGNUM conifer trees and other components of this project have been used for the calculations of publications
- R. Sievänen, J. Perttunen, E. Nikinmaa, and P. Kaitaniemi. Toward extension of a single tree functional- structural model of scots pine to stand level: effect of the canopy of randomly distributed, identical trees on development of tree structure. Functional Plant Biology, 35(9/10):964–975, 2008.
- R. Sievänen, P. Raumonen, J. Perttunen, E. Nikinmaa, and P. Kaitaniemi. A study of crown development mechanisms using a shoot-based tree model and segmented terrestrial laser scanning data. Annals of Botany, 122(3):423–434, 2018.

Other publications that are referred to in this document are
- Perttunen J, Sievänen R, Nikinmaa R, Salminen H, Saarenmaa H, Väkevä J. 1996. LIGNUM: a tree model based on simple structural units. Annals of Botany 77: 87–98.
- Perttunen J, Sievänen R, Nikinmaa E. 1998. LIGNUM: a model combining the structure and the functioning of trees. Ecological Modelling 108: 189–198.


### CMake for macOS and Unix/Linux Makefile build system

To create Makefile build system with CMake first create the
build tree  directory and  then with `cmake`  the Unix  Makefile build
system itself. To build the Lignum core system:

	cd lignum-core
	mkdir build
	cd build 
	cmake .. 
	make install
	
See also *lignum-core* [README](https://github.com/lignumsystem/lignum-core/blob/master/README.md).

To create LignumForest Makefile build system for debug and compile `lignum-forest` binary 
type:

    cd LignumForest
    mkdir debug
    cd  debug
    cmake .. -DCMAKE_BUILD_TYPE=Debug
    make install 

For LignumForest Makefile build system for Release (optimised, no debug information) type:

    cd LignumForest
    mkdir release
    cd release
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make install

In both cases `make install` will move `lignum-forest` to LignumForest directory
where there are two example  shell scripts to run the program:
	
    run-lignum-forest.sh

Command line options and their  short documentation can be obtained by
running `./lignum-forest`  without any  command line parameters.  See also
LignumForest::Usage().

The main growth loop for LignumForest is implemented in *lignum-forest.cc*.

> [!IMPORTANT]
> It is important to type `make install` to also move `lignum-forest` to
> directory above to be used by the scripts to run simulatations.
> Typing just `make` the `lignum-forest` program remains in the compilation directory.


To recompile `lignum-forest` type:

	make clean
	make install
	
> [!IMPORTANT]
> CMake tracks by default file changes only in the current project (e.g. LignumForest in this case). 
> To let CMake  follow all file dependencies correctly `make clean` is mandatory before recompilation. 
> After `make clean` CMake will have correct build tree from previous software  build.

> [!NOTE]
> To remove all CMake  configurations and compilation work just
> remove the build  tree directory (i.e. *debug*,  *release* or *xcode*)
> and recreate the build tree directory.

CMake  projects   are   configured  with   *CMakeLists.txt*
files. For  this CMake  has an  extensive set  of CMake  variables and
built-in functions that can be set in CMakeLists.txt files or given in
command line.

The best way to  learn CMake is by  studying examples.
lignum-core and LignumForest provide  CMakeLists.txt file examples how
to create libraries, find and integrate external libraries (Qt, HDF5),
create and use external binaries (`l2c` to compile L-system files) and
setup the final product with its dependenices.

> [!NOTE]
> It seems Qt4 is becoming difficult maintain in MacPorts. It does not compile on M1 Apple Silicon. 
>Also Qt `qmake` is becoming obsolete; Qt project has switched to CMake since Qt6.

### CMake for Xcode

For Xcode IDE create the Xcode project file:

    mkdir xcode
    cd xcode
    cmake .. -G Xcode

Open  Xcode  IDE  from  Terminal. Alternatively open  the  Xcode  project  file
`lignum-forest.xcodeproj` from XCode:
     
	 open lignum-forest.xcodeproj

Build the `lignum-forest` Product in  Xcode for debugging.  It will appear
in *xcode/Debug*  directory:

	Xcode -> Product (in the menu bar) -> Build For -> Running/Testing/Profiling

See  also that: 

	Xcode -> Product (in the menu bar) -> Scheme 

is set  to `lignum-forest` to allow Run: 

	Xcode -> Product (in the menu bar) -> Run
	
to debug the program. Xcode IDE itself tracks file dependencies.

Copy necessary \*.fun  function files and \*.txt parameter files to
*xcode/Debug*  where   the  `lignum-forest`  is  located   in  this  case.
Otherwise  hard coded  files names  in the  program are  not found  by
`lignum-forest`. You can also copy `lignum-forest` to LignumForest project
directory instead and load the binary to Xcode from there. 

Set command  line parameters for  `lignum-forest` in Xcode:

	Xcode -> Product (in the menu  bar) -> Scheme ->  Edit Scheme -> Arguments.

Divide the command line into practical parts for debugging from `Arguments -> '+'`.

### CMake for LignumForest dependency graph

CMake allows to generate `graphviz` output file to show all library and executable dependencies of the project.
Then with `dot` create image file with desired file format. For example in the *release* directory type:
	
	mkdir graphviz
	cmake ..   --graphviz=graphviz/LignumForest.dot
	dot -Tpdf -Kneato -Goverlap=prism  graphviz/LignumForest.dot  -o  LignumForest.pdf
	
The output file *LignumForest.pdf* contains the visual presentation of the target dependenices including
external binaries and required link libraries. The option `-T` understands many well known image file formats.

### Documentation

The introductionary presentation will appear in [GENERAL_DESCRIPTION](GENERAL_DESCRIPTION.md).

Simulation results are saved in [HDF5 files](HDF5FILES.md).

The Reference Guide for the LignumForest will be based on comments and other information
available in the software. Extraction of the comments, rendition of the software content and 
architecture, generation of the structure of the document and formatting the document to html 
and LaTeX will be done by `doxygen`. To generate the documentation run `doxygen` in LignumForest directory:
    
    doxygen Doxyfile 2> errors.txt
     
Doxyfile is the configuration file for `doxygen`. The documentation will appear in DoxygenDoc directory. 
Errors and warnings will appear in *errors.txt*. To see html version of the document type (on macOS):

    open DoxygenDoc/html/index.html
    
To generate LaTeX version go to latex subdirectory and use make:

    cd DoxygenDoc/latex
    make all
    
The result will be *refman.pdf* that can be opened with a pdf reader.

To use Doxyfile the following three programs are needed:

  + doxygen: generate the document 
  + dot: used by `doxygen` to generate graphs for class hierarchies and function calls.
  + doxywizard: GUI to browse, to edit and optionally to run Doxyfile. 
    
On macOS these are easiest to install with MacPorts.
