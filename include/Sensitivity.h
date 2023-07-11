#ifndef SENSITIVITY_H
#define SENSITIVITY_H
#include <Lignum.h>
///\file Sensitivity.h
///Print tree data for sensitivity analysis
///\deprecated Use HDF5 files for analysis
namespace LignumForest{
  template <class TS,class BUD>
  class Sensitivity{
  public:
    //Close output file 
    ~Sensitivity();
    //Open the output file in append mode and print the header
    void printHeader(const string& file_name);
    //Print data for sensitivity analysis
    void printData(Tree<TS,BUD>& t);
  private:
    ofstream file;
  };
}
#endif

#include <SensitivityI.h>
