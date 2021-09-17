using namespace std;

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <ParseCommandLine.h>
#include <mathsym.h>


void Usage()
{
  cout << "Usage:  ./sens -upfile file -downfile file -nfile file [-scale]" << endl;
}

//Reads three files: nfile, upfile, downfile, in which results of
//simulations have been stored. The parameter values in simulations are
//nominal, higher value (up), lower value (down). Program reads in the
//files and calculates sensitivites of variables (dX/dp) by evaluating
//derivate at nominal parameter value of second order curve going
//through the points (pdown, Xdown), (pnom, Xnom) and (pup,
//Xup). Program reads the name of parameter and value from the first
//line of the file. It is assumed that pdown, pnom and pup are
//different (values). On the second line in the file are names of
//varibles.

//Program writes a file with sensitivities, first line is names of
//variables, ie. column headings. The name of file is formed with the aid of
//the name of parameter.

//If switch -scale is set the individual sensitivitiy is scaled with
//the range of variation (max - min) of the variable.

//Program compiles with command
//> c++ -I../../c++adt/include sens.cc -o senss -L../../c++adt/lib -lcxxadt 

int main(int argc,char* argv[])
{

    if(argc < 5) {
      Usage();
      return -1;
    }

    string upfile;

    if(!ParseCommandLine(argc,argv,"-upfile", upfile ))
      return -1;

    ifstream upf(upfile.c_str());
    if(!upf)
      return -1;

  string line;
  double p_up;
  string par;
  getline(upf,line); //Line with parameter name & value
  istringstream iss(line);
  iss >> par >> p_up; 
  string names; //names of variables
  getline(upf,names);
  double up_v[100][21];
  getline(upf,line);
  int nline = 0;
  while(!upf.eof()){
    istringstream iss(line);
    for(int i=0; i < 21; i++)
      iss >> up_v[nline][i];
    getline(upf,line);
     nline ++;
  }

  upf.close();

  string downfile;

  if(!ParseCommandLine(argc,argv,"-downfile", downfile ))
      return -1;

  ifstream downf(downfile.c_str());
  if(!downf)
    return -1;

  double p_down;
  getline(downf,line); //Line with parameter name & value
  istringstream isd(line);
  isd >> par >> p_down; 
  getline(downf,names);
  double down_v[100][21];
  getline(downf,line);
  nline = 0;
   while(!downf.eof()){
    istringstream iss(line);
    for(int i=0; i < 21; i++)
      iss >> down_v[nline][i];
 
     getline(downf,line);
     nline ++;
   }
   downf.close();

  string nfile;

  if(!ParseCommandLine(argc,argv,"-nfile", nfile ))
      return -1;

  ifstream nf(nfile.c_str());
  if(!nf)
    return -1;

  double p_n;
  getline(nf,line); //Line with parameter name & value
  istringstream isn(line);
  isn >> par >> p_n; 

  getline(nf,names);
  double n_v[100][21];
  getline(nf,line);
  nline = 0;
   while(!nf.eof()){
    istringstream iss(line);
    for(int i=0; i < 21; i++)
      iss >> n_v[nline][i];

 
      getline(nf,line);
     nline++;
   }
   nf.close();

    double s[100][21];

   double x1 = p_down;
   double x2 = p_up;
   double x0 = p_n;

   for(int i = 0; i < nline; i++) {
     for(int j = 1; j < 21; j++) {
       double y1 = down_v[i][j];
       double y2 = up_v[i][j];
       double y0 = n_v[i][j];
       double slope;
       slope = ((x2-x0)*(y1-y0)/(x2-x0)-(x1-x0)*(y2-y0)/(x2-x0))/(x2-x1);
       s[i][j] = slope;
     }
     s[i][0] = n_v[i][0];
   }

   bool scale = false;
   if(CheckCommandLine(argc,argv,"-scale")) 
     scale = true;
   else
     scale = false;

   if(scale) {
     double maxs[21], mins[21];
     for(int j = 1; j < 21; j++) {
       maxs[j] = -R_HUGE;
       mins[j] = R_HUGE;
       for(int i = 0; i < nline; i++) {
	 if(mins[j] > n_v[i][j]) mins[j] = n_v[i][j];
	 if(maxs[j] < n_v[i][j]) maxs[j] = n_v[i][j];
       }
     }

     for(int i = 0; i < nline; i++) {
       for(int j = 1; j < 21; j++) {
	 s[i][j] /= maxs[j]-mins[j];
       }
     }
 
   }

   string out_name = "S_";
   out_name += par;
   out_name += ".dat";
    ofstream out(out_name.c_str(), ofstream::trunc);
   if(out) {
     out << names << endl;
    for(int i = 0; i < nline; i++) {
     for(int j = 0; j < 21; j++) {
       out << s[i][j] << " ";
     }
     out << endl;
    }
   }

   out.close();
} //End of program


