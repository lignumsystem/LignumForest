#ifndef TREELOCATIONS_H
#define TREELOCATIONS_H
#include <iostream>
#include <cmath>
#include <mathsym.h>
#include <vector>
#include <utility> 
///\file TreeLocations.h
using namespace std;
using namespace cxxadt;

namespace LignumForest{
/// Representation of a circular gap: <<x,y>,r>
/// where x,y is the center of gap of raius r
typedef pair<pair<double,double>,double> ForestGap;


void GenerateLocations(int& nTrees, double corner1X, double corner1Y, 
		       double corner2X, double corner2Y, double rmin, 
		       const ForestGap& gap,
		       vector<pair<double,double> >& v);
}
#endif
