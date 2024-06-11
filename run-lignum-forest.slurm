#!/bin/bash
#CHECK --account (Default: $USER) 
#CHECK --mail-user from command line
#SBATCH --job-name=lignum-forest   # create a short name for your job
#SBATCH --account=$USER          # account
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1              # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core
#SBATCH --partition=normal         # sorvi: normal, big or test
#SBATCH --time=90:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --mail-user jari.perttunen@luke.fi
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#
#REMEMBER to update command line, especially the output file name
./lignum-forest -iter 80 -metafile '{MetaFile,MetaFileFip6}.txt' -voxelspace VoxelSpace.txt  -modeChange 10 -architectureChange 10 -numParts 2 -hw 10 -kBorderConifer 0.11 -writeInterval 20 -increaseXi 15 -generateLocations 770 -treeDist 0.3 -hdf5 June11-6-LignumForest.h5  -dumpSelf  -n_buds_ini_min 3 -n_buds_ini_max 3 -verbose
