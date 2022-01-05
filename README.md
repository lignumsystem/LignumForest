# LignumForest
LignumForest is a project for simulating a growing tree community using individual LIGNUM trees.
This project must reside under *lignum-core* directory. That means first
clone [lignum-core](https://github.com/lignumsystem/lignum-core.git) repository and then
in lignum-core clone [LignumForest](https://github.com/lignumsystem/LignumForest.git).

## Compilation
To compile LignumForest (and lignum-core) type:

    cd LignumForest
    qmake  LignumForest.pro
    make

To compile with optimization on (faster, no debug) type:

    qmake  "CONFIG+=release" LignumForest.pro
    make

To remove all compilation work type `make distclean`.
    
## Documentation

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
  + doxywizard: GUI to browse, edit and optionally run Doxyfile. 
    
On macOS these are easiest to install with MacPorts (or some other software package system). 
