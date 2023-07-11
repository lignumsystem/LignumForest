#ifndef MIXED_FOREST_H
#define MIXED_FOREST_H

#include <utility>
#include <tr1/tuple>
#include <GrowthLoop.h>

namespace LignumForest{
  template <class GL1, class GL2>
  class MixedForest{
  public:
    void initialize(int argc, char** argv);
    void growthLoop();
    void afterGrowth();
  private:
    typedef std::tr1::tuple<GL1,GL2> MForest;
    MForest forest;
  };

  template <class GL1,class GL2>
  void MixedForest<GL1,GL2>::initialize(int argc, char** argv)
  {
    for (int i = 0; i < std::tr1::tuple_size<MForest>::value;i++){
      if (i==0)
	std::tr1::get<0>(forest).initialize(argc,argv);
      else if (i==1)
	std::tr1::get<1>(forest).initialize(argc,argv);
    }
  }

  template <class GL1,class GL2>
  void MixedForest<GL1,GL2>::growthLoop()
  {
    int iterations = std::tr1::get<0>(forest).getIterations();
    for(int year = 0; year < iterations; year++) {
      for (int i = 0; i < std::tr1::tuple_size<MForest>::value;i++){
	if (i == 0)
	  std::tr1::get<0>(forest).timeStep(year);
	else if (i == 1)
	  std::tr1::get<1>(forest).timeStep(year);
      }
    }
  }

  template <class GL1,class GL2>
  void MixedForest<GL1,GL2>::afterGrowth()
  {
    for (int i = 0; i < std::tr1::tuple_size<MForest>::value;i++){
      if (i == 0)
	std::tr1::get<0>(forest).afterGrowth();
      else if (i == 1)
	std::tr1::get<1>(forest).afterGrowth();
    }
  }
}//End namespce LignumForest

#endif

