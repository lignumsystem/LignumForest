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

class FindTreeWithProbability{
public:
  bool operator()(const double val,const pair<string,double>& e)const;
};
}
#endif
