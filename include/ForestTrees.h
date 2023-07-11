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
  vector<pair<string,double> > FindTreeFiles(const string& tree_file,double age);
  ///Find a tree with a probability 
  class FindTreeWithProbability{
  public:
    ///\param val Probability to return *true*  for a tree
    ///\param  e Representation of a tree from a file \sa FindTreeFiles
    bool operator()(const double val,const pair<string,double>& e)const;
  };
}
#endif
