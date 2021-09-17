#!/bin/sh
#          iter 
nice -10 ./lig-radiation -iter 20 -metafile MetaFile.txt -voxelspace VoxelSpace.txt -outputAll -Voxbox 0.2 -numParts 5 -treeLocations Radiation/TL-Rad-4.txt -inputTree this-9/this9-1-2-10.xml -phprodfile Valojuttu/ph-prod-1-2-10-Box02.dat -targetTreeRad 0.5 -noWoodVoxel -radMethod 3 > valo-1-2-10-Box02.txt
nice -10 ./lig-radiation -iter 20 -metafile MetaFile.txt -voxelspace VoxelSpace.txt -outputAll -Voxbox 0.2 -numParts 5 -treeLocations Radiation/TL-Rad-4.txt -inputTree this-9/this9-1-2-10.xml -phprodfile Valojuttu/ph-prod-1-2-10-Box01.dat -targetTreeRad 0.5 -noWoodVoxel -radMethod 3 > valo-1-2-10-Box01.txt
nice -10 ./lig-radiation -iter 20 -metafile MetaFile.txt -voxelspace VoxelSpace.txt -outputAll -Voxbox 0.2 -numParts 5 -treeLocations Radiation/TL-Rad-4.txt -inputTree this-9/this9-1-2-10.xml -phprodfile Valojuttu/ph-prod-1-2-10-Box03.dat -targetTreeRad 0.5 -noWoodVoxel -radMethod 3 > valo-1-2-10-Box03.txt
nice -10 ./lig-radiation -iter 20 -metafile MetaFile.txt -voxelspace VoxelSpace.txt -outputAll -Voxbox 0.2 -numParts 5 -treeLocations Radiation/TL-Rad-4.txt -inputTree this-9/this9-1-2-10.xml -phprodfile Valojuttu/ph-prod-1-2-10-Box04.dat -targetTreeRad 0.5 -noWoodVoxel -radMethod 3 > valo-1-2-10-Box04.txt


