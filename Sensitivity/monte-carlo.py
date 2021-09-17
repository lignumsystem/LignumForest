#python script to make 'Monte-Carlo' testing of parameters and functions
#
import sys
import os
import string


#implement the editing of files here

#then run lignum 
result = os.spawnlpe(P_WAIT,"./scotspine","30   2 10 10 MetaFile.txt VoxelSpace.txt -treeDist 0.5 -startVoxCalc 0 -numParts 2 -light 4 -toFile pine-managed.txt -xml Tree30.xml",os.environ)

#if lignum exits with exit(0) the results may be good, move the results aside
if result == 0:
    now = datetime.datetime.now()
    dir = str(now.day)+'-'+str(now.month)+'-'+str(now.hour)+'-'+str(now.minute)
    os.mkdir('LGSim'+dir)
    os.spawnv(P_WAIT,"./move-results",dir)
else:#clean up
    os.remove("SegmentFile.txt")
