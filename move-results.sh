#!/bin/sh
mkdir LGSim/$1
mv  *.xml CrownLimit-*-*-*.txt Locations*.txt  Tree?*.txt VoxelSpace-*.txt  *fipdistrib*.txt LGSim/$1
cp *.fun Tree.txt VoxelSpace.txt LGSim/$1 
cp run-lignum.sh LGSim/$1
