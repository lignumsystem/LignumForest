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
#SBATCH --time=72:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --mail-user jari.perttunen@luke.fi
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#
#REMEMBER to update command line, especially the output file name
#IMPORTANT 1: MetaFiles are used in alphabetical order, use appropriate files names 
#IMPORTANT 2: Make sure the first MetaFile used defines all parameter and function files needed.
#             The second one can define changes only.
#IMPORTANT 3: Protect the Glob expression used for the option -metafile with quotation characters.
#IMPORTANT 4: modeChange is a list of years separated commas,e.g : -modeChange 20,30,40 
#IMPORTANT 4: The option -fsapwdown option is mandatory
#May29-1 with current command line including -segLenVar activating random component. 
./lignum-forest -iter 80 -metafile '{MetaFile,MetaFile1}.txt' -voxelspace VoxelSpace.txt -modeChange 10 -architectureChange 10 -butt_swell_coeff 0.01 -butt_swell_start 0 -segLenVar -numParts 2 -hw 10 -kBorderConifer 0.11 -writeInterval 20 -increaseXi 15 -fsapwdown fsapwdown.fun -generateLocations 770 -treeDist 0.3 -hdf5 Oct20-2-May29-1-buttswell-LignumForest.h5 -dumpSelf -n_buds_ini_min 3 -n_buds_ini_max 3 -verbose 
#./lignum-forest -iter 80 -metafile '{MetaFile,MetaFile1,MetaFileLGPmin}.txt' -voxelspace VoxelSpace.txt -modeChange 20,60 -architectureChange 20 -numParts 2 -hw 10 -kBorderConifer 0.11 -writeInterval 20  -fsapwdown fsapwdown.fun -increaseXi 15 -generateLocations 770 -treeDist 0.3 -hdf5 Sep30-4-May29-1-LignumForest.h5 -dumpSelf -n_buds_ini_min 3 -n_buds_ini_max 3 -verbose 
#./lignum-forest -iter 80 -metafile '{MetaFile,MetaFileFip6}.txt' -fsapwdown fsapwdown5.fun -voxelspace VoxelSpace.txt -modeChange 10 -architectureChange 10 -numParts 2 -hw 10 -kBorderConifer 0.11 -writeInterval 20 -increaseXi 15 -xiIncrement 15.0 -generateLocations 770 -treeDist 0.3 -hdf5 August28-5-LignumForest.h5 -dumpSelf -n_buds_ini_min 3 -n_buds_ini_max 3 -verbose
#./lignum-forest -iter 80 -metafile '{MetaFile,MetaFileFip6}.txt' -fsapwdown fsapwdown3.fun -voxelspace VoxelSpace.txt -modeChange 10,20 -architectureChange 10 -Lmaxturn 80 -numParts 2 -hw 10 -kBorderConifer 0.11 -writeInterval 20 -increaseXi 15 -xiIncrement 15.0 -generateLocations 400 -treeDist 0.3 -hdf5 August29-5-LignumForest.h5 -dumpSelf -n_buds_ini_min 3 -n_buds_ini_max 3 -verbose


