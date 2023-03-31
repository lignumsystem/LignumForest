#ifndef FINDFORESTBOUNDINGBOX_H
#define FINDFORESTBOUNDINGBOX_H
#include <Lignum.h>
using namespace std;
namespace LignumForest{ 
template <class TREE, class TS, class BUD>
class ForestBoundingBoxIndividualTrees{
public:
  BoundingBox operator()(BoundingBox bb, TREE* t)
  {
    //Find the bounding box for a single tree
    BoundingBox base;
    base = Accumulate(*t,base,FindCfBoundingBox<TS,BUD>());
    //Adjust the two opposite corners of the bounding box if necessary
    //Move the 'origo', lower left hand corner
    bb.setMinX(min(bb.getMin().getX(),base.getMin().getX()));
    bb.setMinY(min(bb.getMin().getY(),base.getMin().getY()));
    bb.setMinZ(min(bb.getMin().getZ(),base.getMin().getZ()));
    //Move the opposite corner, the upper right hand corner
    bb.setMaxX(max(bb.getMax().getX(),base.getMax().getX()));
    bb.setMaxY(max(bb.getMax().getY(),base.getMax().getY()));
    bb.setMaxZ(max(bb.getMax().getZ(),base.getMax().getZ()));
    return bb;
  }
};

//Find the bounding box for a single tree 
template <class TS, class BUD>
class ForestBoundingBoxTreeCopies{
public:
  ForestBoundingBoxTreeCopies(Tree<TS,BUD>& t):tree(t){}
  BoundingBox operator()(BoundingBox bb,pair<double,double> p1)
  {
    //Move the tree to next location
    const Point& p2 = GetPoint(tree);
    MoveTree<TS, BUD> move(Point(p1.first-p2.getX(),
				 p1.second-p2.getY(),0.0),
			   tree);
    ForEach(tree, move);

    //Find the bounding box for a single tree
    BoundingBox base;
    base = Accumulate(tree,base,FindCfBoundingBox<TS,BUD>());
    //Adjust the two opposite corners of the bounding box if necessary
    //Move the 'origo', lower left hand corner
    bb.setMinX(min(bb.getMin().getX(),base.getMin().getX()));
    bb.setMinY(min(bb.getMin().getY(),base.getMin().getY()));
    bb.setMinZ(min(bb.getMin().getZ(),base.getMin().getZ()));
    //Move the opposite corner, the upper right hand corner
    bb.setMaxX(max(bb.getMax().getX(),base.getMax().getX()));
    bb.setMaxY(max(bb.getMax().getY(),base.getMax().getY()));
    bb.setMaxZ(max(bb.getMax().getZ(),base.getMax().getZ()));
    return bb;
  }
private:
  Tree<TS,BUD>& tree;
};


//Resize the bounding box 
class ResizeBoundingBox{
public:
  ResizeBoundingBox(const Point& tp, BoundingBox& b):tree_point(tp),box_tree(b){}
  BoundingBox operator()(BoundingBox bb,pair<double,double> p)
  {
    Point minimum = box_tree.getMin();
    Point maximum = box_tree.getMax();

    //Resolve the lower  left corner
    double add_x = p.first-tree_point.getX();
    double add_y = p.second-tree_point.getY();
    //Move the lower left corner
    Point p1(minimum.getX()+add_x,minimum.getY()+add_y,minimum.getZ());
    //Move the upper right hand corner
    Point p2(maximum.getX()+add_x,maximum.getY()+add_y,maximum.getZ());

//     cout << bb<<flush;
//     cout << box_tree <<flush;
//     cout << "Tree " << tree_point<<flush;
//     cout << "add_x " << add_x << " add_y " << add_y <<endl;
//     cout << "P    " << p.first << " "<< p.second <<endl;
//     cout << "p1   " << p1 << flush;
//     cout << "p2   " << p2 <<flush;
    
    //Adjust the two opposite corners of the bounding box if necessary
    //Move the 'origo', lower left hand corner
    bb.setMinX(min(bb.getMin().getX(),p1.getX()));
    bb.setMinY(min(bb.getMin().getY(),p1.getY()));
    //Move the opposite corner, the upper right hand corner
    bb.setMaxX(max(bb.getMax().getX(),p2.getX()));
    bb.setMaxY(max(bb.getMax().getY(),p2.getY()));
//     cout << bb<<endl;
    return bb;
  }
  const Point& tree_point; 
  BoundingBox& box_tree;
};
//Find the bounding box  for conifer forest applying FindCfBoundingBox
//for a single tree and then using the bounding box found and the tree
//locations resize the found bounding box 
template <class TS, class BUD>
BoundingBox FindForestBoundingBox(Tree<TS,BUD>& tree, const vector<pair<double,double> >& locations)
{
  BoundingBox box;
  box = Accumulate(tree,box,FindCfBoundingBox<TS,BUD>());
  box = accumulate(locations.begin(),locations.end(),box,ResizeBoundingBox(GetPoint(tree),box));
  return box;
}

//Find the bounding box for conifer forest applying FindCfBoundingBox for a single tree.
template <class TREE, class TS, class BUD>
BoundingBox FindForestBoundingBox(vector<TREE*>& vtree, vector<pair<double,double> >& locations,
				  const Point& lower_left, const Point& upper_right)
{
  BoundingBox bb;
  //Many individual trees 
  if (vtree.size() > 1){
    //Initialize lower left and upper right hand corners
    bb.setMin(lower_left);
    bb.setMax(upper_right);
    //Find the bounding box for each tree and adjust forest bounding box accordingly
    bb = accumulate(vtree.begin(),vtree.end(),bb,ForestBoundingBoxIndividualTrees<TREE,TS,BUD>());
  }
  //One single tree in the forest
  else if (vtree.size() == 1){
    //Save the original position of the tree
    Point p1 = GetPoint(*vtree[0]);
    //Initialize lower left and upper right hand corners
    bb.setMin(lower_left);
    bb.setMax(upper_right);
    //Find the forest bounding box by moving the target tree from one position to another
    bb = accumulate(locations.begin(),locations.end(),bb,ForestBoundingBoxTreeCopies<TS,BUD>(*vtree[0]));
    Point p2 =  GetPoint(*vtree[0]);
    //Move the tree back to its original position
    MoveTree<TS, BUD> move(p1-p2,*vtree[0]);
    ForEach(*vtree[0], move);
  }
  else{
    cout << "FindForestBoundingBox error in number of trees " << vtree.size() <<endl;
  }
  return bb;
}
}
#endif
