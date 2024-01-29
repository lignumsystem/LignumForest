# Tasks before simulating Scots pine with LignumForest

CrownDensity simulates Scots pine stand in satisfactorily. 
Before simulating Scots pine stand with LignumForest
some program changes and checks are needed.

## Task list
- [ ] Change LignumForest command line identical to CrownDensity command line.
- [ ] Implement required changes by the new command line in the main loop
- [x] Check L-system is identical to L-system in CrownDensity
- [ ] Change SetSegmentLength to one used in CrownDensity
- [ ] Check the Meta files, parameter and function files in CrownDensity
      that they the ones that produce a satisfactory tree development
- [ ] Change Meta files, parameter and function files to ones 
      used in CrownDensity
- [ ] Transfer global variables used from CrownDensity to LignumForest
      and use LignumForest namespace
	  
## Additional tasks

It seems CGAL 5.6 (latest manual in cgal.org) can write 
Paraview *vtp* files and also *.obj* used by Blender. 
This could help transition from LignumWb to Paraview and/or Blender 
or perhaps even relatively straightforward at best.

- [ ] Study how to create segment cylinders and leaf polygons in CGAL.
      These are probably repreented as polygon meshes or some other 
	  low level geometric primitives supported by CGAL.

