#ifndef TREELOCATIONS_H
#define TREELOCATIONS_H
#include <iostream>
#include <cmath>
#include <mathsym.h>
#include <vector>
#include <utility> 
///\file TreeLocations.h
///\brief Genrate forest plot tree locations
using namespace std;
using namespace cxxadt;

namespace LignumForest{
  /// Representation of a circular gap: <<x,y>,r>
  /// where x,y is the center of gap of raius r
  typedef pair<pair<double,double>,double> ForestGap;

  /// \brief Generate tree locations
  /// Generate tree loacations as Nonstationary Poisson process: Stoyan, Kendall & Mecke, p. 52 -54
  /// \param[in,out] nTrees Number of trees as target number in a forest plot and number of trees actually generated.
  /// \param corner1X Lower left corner X
  /// \param corner1Y Lower left corner Y
  /// \param corner2X Upper right corner X
  /// \param corner2Y Upper right corner Y
  /// \param rmin     Minimum distance between trees
  /// \param gap      Representation of a circular gap. ForestGap is simply pair<pair<x,y>,r> 
  ///                 where pair<x,y> is the position and r the radius of the circular gap. \sa ForestGap
  /// \param[in,out] v       Vector to contain generated tree locations
  void GenerateLocations(int& nTrees, double corner1X, double corner1Y, 
			 double corner2X, double corner2Y, double rmin, 
			 const ForestGap& gap,
			 vector<pair<double,double> >& v);
}

#endif
