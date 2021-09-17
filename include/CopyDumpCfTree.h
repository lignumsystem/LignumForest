#ifndef COPYDUMPCFTREE_H
#define COPYDUMPCFTREE_H

#include <Lignum.h>
#include <VoxelSpace.h>

//For a sequence (vector) of tree  locations move and dump the tree to
//specified locations in the voxel space.
template <class TS, class BUD>
class CopyDumpScotsPineTree{
public:
  CopyDumpScotsPineTree(Tree<TS,BUD>& t, VoxelSpace& vs, int n):tree(t),voxel_space(vs),num_parts(n){}
  void operator()(pair<double,double>& p)const{
    //Move the tree to next location
    const Point& point = GetPoint(tree);
    MoveTree<TS, BUD> move(Point(p.first-point.getX(),
				 p.second-point.getY(),0.0),
			   tree);
    ForEach(tree, move);
    //dump the Scots pine tree into the voxelspace
    DumpScotsPineTree(voxel_space,tree,num_parts);
  }
private:
  Tree<TS,BUD>& tree;
  VoxelSpace& voxel_space;
  int num_parts;
};

#endif
