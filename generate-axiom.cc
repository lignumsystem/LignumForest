#include <GenerateAxiom.h>
#include <XMLTree.h>
#include <Lignum.h>
#include <ScotsPine.h> 
#include <cstring>
#include <fstream>

using namespace std;

//Includes all kinds of stuff, turtle graphics etc.
#include <lengine.h>

int ran3_seed;

double H_0_ini, H_var_ini;          //For variation of initial heights and
int n_buds_ini_min, n_buds_ini_max;   // and number of buds (in .L file)
double  rel_bud;
bool  bud_variation;
double branch_angle;

//and for pine
namespace Pine{
#include <LSystem.h>

}



void Usage()
{
  cout << "Usage: ./generate-axiom -xml <input_file> -axiom <output_file> [-test]" 
       << endl;
  exit(0);
}

//A  program to  generate Axiom  from the  existing tree.  Include the
//output file to your existing L system as an axiom. Usage:

//   ./generate-axiom -xml <input_file> -axiom <output_file> [-test]

//At this moment also three  output tree files (XML) are generated: 1)
//the original  tree file, 2)  tree after lstringToLignum and  3) tree
//after  lignumToLstring.  The  axiom generation  is designed  for the
//Scots  pine used  in ForestLignum.   The testing  includes  only the
//axiom  itself (l.start()).   But the  subsequent  L-string expansion
//should work as usual.

//The project file is  GenerateAxiom.pro and includes the pine-axiom.L
//as a  test case.  You  may want test  your own axiom instead  of the
//current  one.  To do  this you  have to  compile the  project twice,
//starting from the  qmake command. First, the new  axiom is generated
//(without  -test  option)  that   has  to  be  manually  inserted  to
//pine-axiom.L . Secondly, this new L-system is compiled and linked to
//this program and you can test the axiom with the -test option.
int main(int argc,char* argv[])
{
  string input_file_name;
  string output_file_name;
  fstream ifile;
  fstream ofile;
  
  ParseCommandLine(argc,argv,"-xml",input_file_name);
  if (input_file_name.empty())
    Usage();
  
  ParseCommandLine(argc,argv,"-axiom",output_file_name);
  if (output_file_name.empty())
    Usage();

  ifile.open(input_file_name.c_str(),ios_base::in);
  ofile.open(output_file_name.c_str(),ios_base::out);

  Tree<ScotsPineSegment,ScotsPineBud> pine(Point(0,0,0), 
					   PositionVector(0,0,1));
  Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData> l;
  XMLDomTreeReader<ScotsPineSegment, ScotsPineBud> reader;
  XMLDomTreeWriter<ScotsPineSegment, ScotsPineBud> writer;

  reader.readXMLToTree(pine,input_file_name);

  //GenerateAxiom generates the axiom that can be included to an existing LSystem
  cout << "Generating the axiom " << endl;
  GenerateAxiom(ofile,pine);
  cout << "Done. NOTE: Please include declarations for" << endl
       << "SetHead(double,double,double) and MoveTo(double,double,double)" << endl
       << "modules in your L-system file"<< endl;
  
  //You can remove/uncomment this part if you will. But this is also a
  //handy way to test the axiom
  if (CheckCommandLine(argc,argv,"-test")){
    cout << "Testing and generating three files" << endl;
    cout << "1. The original tree to file " << "axiom-pine-1.xml" <<endl;
    writer.writeTreeToXML(pine,"axiom-pine-1.xml");
    cout << "2. Expanding axiom" <<endl;
    l.start();
    cout << "3. lstringToLignum to file " << "axiom-pine-2.xml" <<endl;
    l.lstringToLignum(pine,1,PBDATA);
    writer.writeTreeToXML(pine,"axiom-pine-2.xml");
    cout << "4. lignumToLstring to file "  << "axiom-pine-3.xml" << endl;
    l.lignumToLstring(pine,1,PBDATA);
    writer.writeTreeToXML(pine,"axiom-pine-3.xml");
    cout << "Done" << endl;
  }
  exit(0);
}
