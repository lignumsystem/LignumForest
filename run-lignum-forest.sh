#CHECK: Have you done `make install` after compilation?
#REMEMBER: Set -hdf5 <filename> !!!!!!!
#NOTE: Protect regular expression for -metafile with quotes. 
./lignum-forest -iter 22 -metafile MetaFile.txt -voxelspace VoxelSpace.txt -butt_swell_coeff 0.0015 -butt_swell_start 0 -numParts 2 -hw 5 -kBorderConifer 0.14 -writeInterval 20 -increaseXi 15 -fsapwdown fsapwdown.fun -generateLocations 770 -treeDist 0.3 -hdf5 koe1.h5 -dumpSelf -n_buds_ini_min 3 -n_buds_ini_max 3 -verbose -analyze_k 20 
