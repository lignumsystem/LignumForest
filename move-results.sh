#!/bin/sh
mkdir LGSim/$1
mv  *.h5 *.pdf   LGSim/$1
scp *.fun Tree*.txt VoxelSpace.txt LGSim/$1 
scp run-lignum-forest.slurm LGSim/$1
