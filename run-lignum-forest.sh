#CHECK: Have you done `make install` after compilation?
#REMEMBER: Set -hdf5 <filename> !!!!!!!
#NOTE: Protect regular expression for -metafile with quotes. 
./lignum-forest -iter 60 -metafile 'MetaFile*.txt' -voxelspace VoxelSpace40x40.txt  -modeChange 20 -architectureChange 10 -numParts 2 -hw 10 -kBorderConifer 0.11 -writeInterval 5 -increaseXi 15 -generateLocations 4000 -treeDist 0.3 -hdf5 May162024LignumForest.h5  -dumpSelf  -n_buds_ini_min 3 -n_buds_ini_max 3  -verbose > runlog.txt
