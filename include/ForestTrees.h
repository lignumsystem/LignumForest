/// \file ForestTrees.h
/// \brief Find tree files 
#ifndef FORESTTREES_H
#define FORESTTREES_H

#include <QDir>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <utility>
using namespace std;
namespace LignumForest{
  ///\brief Find tree files
  ///
  ///Given the  \p tree_file parse it and  return file names  that match and
  ///exist in the current directory, associate the given probability with
  ///the file name
  ///\param tree_file Tree file names
  ///\param age Tree age
  ///\return Vector of tree files with associated probability
  vector<pair<string,double> > FindTreeFiles(const string& tree_file,double age);
  ///\brief Find a tree with a probability 
  class FindTreeWithProbability{
  public:
    ///\param val Probability to return *true*  for a tree
    ///\param  e Representation of a tree from a file
    ///\sa FindTreeFiles
    bool operator()(const double val,const pair<string,double>& e)const;
  };
}
#endif
