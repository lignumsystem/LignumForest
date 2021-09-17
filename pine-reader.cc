//Include Lignum implementation 
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <Lignum.h>
#include <Bisection.h>
#include <Shading.h>
#include <TreeLocations.h>
#include <CopyDumpCfTree.h>
//XML file 
#include <XMLTree.h>
//Include the implementation of the tree segment and bud
#include <ScotsPine.h>
#include <CalculateLight.h>
#include <VoxelSpace.h>
#include <VisualFunctor.h>

//Impelements VisualizeLGMTree
#include <GLSettings.h>
#include <OpenGLUnix.h>
#include <LGMVisualization.h>

//Includes all kinds of stuff, turtle graphics etc.
#include <lengine.h>

//and for pine, see also pine9bp.L in lsys.
namespace Pine{
#include <LSystem.h>

}

#include <HarvestStand.h> 
#include <SomeFunctors.h>
#include <DiameterGrowth.h>

void Usage()
{
  cout << "Usage:  ./scotspine-reader -xml <filename> [-viz]" <<endl;
  exit(0);
}

int main(int argc, char** argv)
{
  ScotsPineTree pine1(Point(0,0,0),PositionVector(0,0,1),
		      "sf.fun","fapical.fun","fsapwdown.fun");

  string xmlfile;
  string clarg;
  if (ParseCommandLine(argc,argv,"-xml",clarg))
    if (clarg.length() > 0)
      xmlfile = clarg;
    else
      Usage();

  //Read the tree from xml file
  if (!xmlfile.empty()){
    XMLDomTreeReader<ScotsPineSegment,ScotsPineBud> writer;
    writer.readXMLToTree(pine1,xmlfile);
  }

  if(CheckCommandLine(argc,argv,"-viz")) {
    LGMVisualization viz;
    viz.InitVisualization(argc,argv);
    // textures 512x512
    viz.AddCfTree(pine1, "Manty.bmp", "neulaset5.tga");
    float th = (float)GetValue(pine1,LGAH);
    cout << th << endl;
    //viz.ResetCameraPosition(th);
    viz.SetMode(SOLID);
    //viz.ResetCameraPosition(GetValue(pine1,LGAH));
    viz.StartVisualization();
  }
 
}
