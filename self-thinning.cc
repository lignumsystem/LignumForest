#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;
//Self thinning model, d > 12cm
double fN12(double d)
{
  if (d == 0.0)
    //Special case d = 0.0
    return fN12(12.0) + 12.0*1.127*pow(10.0,3.0);
  else
    return 1.08685*pow(10.0,6.0)/pow(d,2.0551);
}

//Input 
//Dbase: diameter at base (in meters)
//Dbh:   diameter at 1.3 (in meters)
//N: trees per ha
//Output: trees per ha after self thinning
double SelfThinning(double Dbase, double Dbh, double N)
{
  double deltaN = 0.0;
  //meters to centimeters
  double d = max(2.0+1.25*Dbh*100.0,Dbase*100.0);
  if (d < 12.0){
    //thinning according to the line  dfN/dD (at point D 12cm dfN/dD =
    //-11510.0)
    if (N < fN12(0.0) - 1127.0*d)
      deltaN = 0.0;
    else
      deltaN = N - fN12(0.0)+1127.0*d;
  }
  else{
    //self thinning according to fN
    if (N < fN12(d))
      deltaN = 0.0;
    else
      deltaN = N - fN12(d);
  }

  return N - deltaN;
}

#ifdef SELFTHINNING
int main()
{
  double N = fN12(0.0);
  //  cout <<  0.0 << " "  << N <<endl;
  for (double d = 0.01; d < 0.40;){
    N = SelfThinning(d,d,N);
    //    cout << d << " "  << N << endl;
    d = d + 0.01;
  }
  return 0;
}
#endif
