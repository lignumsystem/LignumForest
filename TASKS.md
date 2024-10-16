# Tasks before simulating Scots pine with LignumForest

CrownDensity simulates Scots pine stand in satisfactorily. 
Before simulating Scots pine stand with LignumForest
some program changes and checks are needed.

## Task list
> [!NOTE]
> It seems all items below are now checked. LignumForest can be tested in sorvi server.

### Command line

	- [x] Added script run-lignum-forest.sh to be used for lignum-forest binary.
	- [x] Merge CrownDensity command line with LignumForest command line.
	- [x] Command line of the run-lignum-forest.sh copied from run-crowndens-basic-model.sh.
	- [x] LignumForest specific arguments present (for example -generateTrees, -treeDist)

### LignumForest main growth loop 

	- [x] Risto/Jari: Check the main growth loop (lignu-forest.cc and \ref AMAIN) is what is required, especially
	- [x] LignumForest::GrowthLoop::createNewSegments() is what is required. It has more that just creating segments. OK?
	- [x] LignumForest::GrowthLoop::allocationAndGrowth is what is required, It has (behing boolean flags) more than pipe model. OK?
	- [x] Jari: Check LignumForest::GrowthLoop::initializeTrees() is what is required, especially
	- [x] LignumForest::branch_angle set to 45 degrees. Branch angle also set in L-system. OK?
	- [x] LignumForest::bud_view_f set also here

### Synchronizing with CrownDensity 

	- [x] Implementation required changes by the new command line in the main loop.
	  
Check and implement the following command line options: <br>
	
	- [x] Jari: Implement checking and setting growth mode change and checking and triggering  architecture change year in the main loop.
	- [x] "-iter"
	- [x] "-metafile"
	- [x] "-voxelspace"
	- [x] "-voxelCalculation"  No "-voxelCalculation" in `lignum-forest`, removed from command line. VoxelSpace always in use. See "-pairwiseSelf"*. 
	- [x] -modeChange
	    - [x] Wild card search for MetaFile*.txt
		- [x] GrowthLoop::growthModeChange in the main loop
	- [x] -architectureChange  *Command line argument implemented, global variables triggering L-system architecure change set*. 
	- [x] -numParts
	- [x] -hw 
	- [x] -kBorderConifer
	- [x] -writeInterval
	- [x] -increaseXi *The same implementation as in CrownDensity*
	- [x] -generateLocations
	- [x] -treeDist 
	- [x] -hdf5
	- [x] Check L-system is identical to L-system in CrownDensity
	- [x] Change SetSegmentLength to one used in CrownDensity. *Using explicitely LignumForest::SetScotsPineSegmentLength*
	- [x] Check the Meta files, parameter and function files in CrownDensity that they are the ones that produce a satisfactory tree development
	- [x] Change Meta files, parameter and function files to ones  used in CrownDensity. Especially -modeChange requires rereading MetaFiles.
	- [x] Transfer global variables used from CrownDensity to LignumForest and use LignumForest namespace.  *LignumForest compiles*. 

### Additional tasks

It seems CGAL 5.6 (latest manual in cgal.org) can write 
Paraview VTK (VTU / VTP) files and also Wavefront Advanced Visualizer 
Object Format (OBJ) used by Blender. This could help transition from 
LignumWb to Paraview and/or Blender or perhaps even relatively 
straightforward at best.

	- [ ] Study how to create segment cylinders and leaf polygons in CGAL.
      These are probably represented as polygon meshes or some other 
	  low level geometric primitives supported by CGAL.

The manual page [CGAL 5.6 - IO Streams]( https://doc.cgal.org/latest/Stream_support/index.html)
might be a good starting point.
