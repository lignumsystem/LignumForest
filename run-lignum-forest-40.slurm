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
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --mail-user jari.perttunen@luke.fi
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
./lig-forest -iter 40 -metafile MetaFile.txt -voxelspace VoxelSpace.txt  -writeOutput -outputAll -writeInterval 1 -numParts 2  -toFile forest40/OK.dat -generateLocations 10 -dumpSelf -hw 15 -increaseXi 15 -xml forest40/eka.xml -verbose -EBH -budViewFunction -adHoc
