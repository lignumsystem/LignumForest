# General Program Description

To appear. This file contains general description of the program and (automatically created) 
links to required libraries and their documentations in The LIGNUM System.
<br> <br>

Run the program without command line options <CODE> ./lig-forest </CODE> to see this

    ./lig-forest -iter <value>  -metafile <file>  -voxelspace <file>
    [-numParts <parts>]  [-treeDist <dist>] [-hw <hw_start>] [-viz]
    [-toFile <filename>] [-xml <filename>] [-writeVoxels] [-sensitivity <filename>]
    [-fipdistrib <filename>] [-writeInterval interval]
    [-seed <num>] [-increaseXi <value>] [-targetTree <num>]
    [-treeFile <filename>] [-generateLocations  <num>] [-woodVoxel] [-treeLocations <file>]
    [-writeOutput] [-verbose] [-bracketVerbose] [-noBorderForest] [-seed <value>] [-kBorderConifer <value>]
    [-H_0_ini <value>] [-H_var_ini <value>] [-n_buds_ini_min <num>] [-n_buds_ini_max <vlaue>]
    [-p0Var <value>] [-segLenVar <value>] [-pairwiseSelf] [-budVariation <value>] [-eero]
    [-gFunVar <value>] [-branchAngleVar <value>]
    [-space0] [-space1] [-space2] [-adHoc]
    [-budViewFunction] [-EBH] -EBH1 <value>]
    [-space2Distance <Value>]
    
    -generateLocations <num>  In this case <num> trees will be generated to random locations. If this
              is not on, tree locations will be read from file Treelocations.txt. This file can be changed
              by -treeLocations <file>. If location file is not found program stops.
    -woodVoxel                If woody parts are dumped to voxels (default = true)
    -treeDist <dist>          Minimum distance between two trees (default = 0), works only with -generateLocations.  
    -numParts <parts>         Segments can be dumped to voxels by parts (i.e. they may belong to different voxels,
                              default = 1)
    -hw <hw_start>            Starting year of formation of sapwood. Default = 15 years
    -increaseXi               If parameter ksi (fraction of heartwood in new segments) increases, value = starting year (default = 15)
                              as increase = 0.004/year after year 15 up to value 0.85 (as in FPB 35: 964-975 2008 article)"
    -targetTree <num>         Any one of the trees can be identified as target tree (default = 0)
    -writeOutput              Most of the things are written to their respctive file at -writeInterval interval (default false)
      -verbose                Output of progress of run if set (default = not set).
    -bracketVerbose           If set, iteration information is printed out in allocation (default = not set).
    -noBorderForest           No border forest around the stand (default = there is border forest)
    -seed <value>             seed for random number generator.
    -kBorderConifer <value>   Extinction coefficient for conifer foliage in border forest (default = 0.14)
    -H_0_ini, -H_var_ini      For variation of initial heights (defaults = 0.3 and 0.0)
    -n_buds_ini_min, -n_buds_ini_max  For variation of initial number of buds (defaults = 4 and 4)
    -p0Var <value>            Random variation of p0 +- max <value> per cent from the value in Tree.txt
    -segLenVar <value>         Random variation of length of new segments around Lnew, per cent
    -pairwiseSelf      Pairwise radiation calculation for the tree itself.
    -eero              For studying the relationships between a) leaf light climate and b) syncronized variation
                       in leaf nitrogen concentration, leaf mass per area and leaf longevity, shoot length and
                       shoot leaf area and c) axis thickness scaling from the base of the stem to the axis tips.
    EBH resource distn can be in use in two ways. Both are set by command line arguments.
    -EBH               means EBH is in use and values (of lambda parameter) are
                       specified for all Gravelius orders (function SPEBHF, ScotsPine.h, function
                       file is specified the constructor of the tree).
    -EBH1 <value>      means that EBH is in use and one value <value> is used
                       for all Gravelius orders. Option -EBH1 <value> overrides option -EBH
                       (EBH is set as SPis_EBH (Scots Pine Parameter Double SPPD), thus 0 == false, 1 == true)
                       EBH is according to W. Palubicki and K. Horel and S. Longay and
                       A. Runions and B. Lane and R. Mech and P. Prusinkiewicz. 2009.
                       Self-organizing tree models for image synthesis ACM Transactions on Graphics 28 58:1-10.


