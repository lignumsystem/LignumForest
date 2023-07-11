#include <iostream>
using namespace std;
//Defines ran3 etc
#include <mathsym.h>
using namespace cxxadt;
#include <LGMUnits.h>
using namespace Lignum;
#include <Pine.h>
using namespace PineTree;
//Defines headers necessary for lsystem
#include <lengine.h>

#include <CrownDensityGlobals.h>
using namespace CrownDensity;



//Start the namespace (compiles to "namespace Pine{")
//The matching "close" statement closes the namespace 
open Pine;

//Defines the lsystem rewriting process. 
#include <lsys.h>
//Global variables
//double L_age = L_H = 0.0;
//a guess for the segment shortening
//const double l1 = 0.9;


const double PI = 3.1415926535897932384;
const double pitch = 45.0*PI/180.0; //the angle of branching in the main axis    = 45
const double roll  = 137.5*PI/180.0; //137.5 is convenient to distribute buds (137.5*N never overlaps)
const double turn_branch = 30.0*PI/180.0;//the angle of branching in branches  = 30
const double pitch_down = 2.0*PI/180.0;//branches bending down annually
const double max_pitch_down = 90.0*PI/180.0;//horizontal branches
const double hupdir = 85.0*PI/180.0; //if a bud points down make it N degrees upwards from horizontal
const double hdir   =  45.0*PI/180.0; //branch buds maintain N degree angle  up from horizontal

const double rollmain2  =(360.0/2.0)*PI/180.0;
const double rollmain3  =(360.0/3.0)*PI/180.0;
const double rollmain4  =(360.0/4.0)*PI/180.0;
const double rollmain5  =(360.0/5.0)*PI/180.0;
const double rollmain6  =(360.0/6.0)*PI/180.0;
const double rollmain7  =(360.0/7.0)*PI/180.0;
const double rollmain8  =(360.0/8.0)*PI/180.0;

//L system runs in two modes. first create the new segments and
//after their lengths are adjusted, the second mode creates the new buds.
//There are two calls to derive in the main program in pine.cc
//NOTE for pine-forest: As the 'mode' is a global variable, it is visible to a l l L-systems.
//Consequently each l[i].derive() will change the mode and alternately the L-systems 
//will receive 0 and 1 as the mode. This is not what we want but what is needed is an additional
//explicit call to l[i].endEach() after each l[i].derive() in pine-forest.cc that resets the mode.
int mode = 0; //0 = new segment, 1 = new buds

//Number of buds as the function of foliage mass
ParametricCurve fnbuds("fnbuds.fun");

//Adjust the number of buds, more light more buds, less light less buds
ParametricCurve fnbudslight("fnbudslight.fun");
//Fipi for buds, the death zone
ParametricCurve fipbud("fip-bud.fun");
module F(double);
//With this implementation side branches start physiological aging, i.e. increments of physiological age.
//When physiological age reaches mode_change_year the side branches create only two new side branches
//Main axis retain physiological age 1
//Bud(data, omega, initial length,physiological age)
module B(PineBudData,int,double,double);
module Pitch(double);
module Roll(double);
module Turn(double);
module HDir(double);
module HUp(double);
module Down(double); 

derivation length: 9;

Start:
{
// If L-system is initilized several times in one program, the functions
// are read in with every initilization
   fnbuds.install("fnbuds.fun");
   fnbudslight.install("fnbudslight.fun"); 
   fipbud.install("fip-bud.fun");

  PineBudData data(ALIVE,0.0,1.0,1.0,1.0);
  mode = 0;
  produce F(0.30) SB() Down(1.0) Pitch(pitch) B(data,2,l1,1.0) EB() 
                  SB() Roll(rollmain4) Down(1.0) Pitch(pitch)  B(data,2,1.0,1.0) EB()
	          SB() Roll(2.0*rollmain4) Down(1.0) Pitch(pitch) B(data,2,1.0,1.0) EB()
		  SB() Roll(3.0*rollmain4) Down(1.0) Pitch(pitch)  B(data,2,1.0,1.0) EB()
                  B(data,1,1.0,1.0);
}

StartEach:
{
}

EndEach:
{
  mode = (mode+1) % 2;
} 

B(data,omega,l,phys_age):
{
  if (mode == 0){
   if (data.state == ALIVE){
     //Note there are no start and end branches
     //These are genrerated in mode 1
     if (omega==1)
       produce F(l) Roll(roll) B(data,omega,l,phys_age);
     else
       produce F(l) B(data,omega,l,phys_age);
   }
   else{
       produce B(data,omega,l,phys_age);
   }
  }
  else{//mode == 1
    //Tassa mukana myos kokeilu tiheydella latvuksen sisalla
    //tiheyden funktio modifioi silmujen lukumaaraa
	
    //first, use the fnbuds for the number of buds
    int nbuds = static_cast<int>(fnbuds(data.fm));
    //then adjust the number of buds with relative light  
    //nbuds = static_cast<int>(max(nbuds+fnbudslight(data.ip),0.0));
    //cout << "ip " << data.ip << " fip(ip) " << fip(data.ip) 
    //     << " FM " << data.fm << " fnbuds.fun " << nbuds <<endl; 

    if(fipbud(data.ip) == 0.0) nbuds = 0;
//    if (nbuds == 1 || nbuds == 2){
//      if (fipbud(data.ip) < 0.6){
//        nbuds = 0;
//      }
//    }

//    tiheysmodifikaatio
      if (is_bud_view_function) {
         nbuds = static_cast<int>(static_cast<double>(nbuds) * bud_view_f(data.view));
      }

    if (nbuds == 0 && L_age < 4)
       nbuds = 1;

    if (data.state == DEAD){
       //cout << "Already dead " << endl;
       produce B(data,omega,l,phys_age);
    }
    else if (nbuds == 0){
       data.state = DEAD;
       //cout << "No new buds, DEAD" <<endl;
       produce B(data,omega,l,phys_age);
    }
    else if (omega == 1){//Main_axis, do not increase physiological age
       if (nbuds == 1)
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1,phys_age) EB()
                   B(data,omega,l,phys_age);
       else if (nbuds == 2)
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1,phys_age+1) EB() 
                   SB() Roll(rollmain2) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
                   B(data,omega,l,phys_age);
       else if (nbuds == 3)
          produce  SB() Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB() 
                   SB() Roll(rollmain3) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
	           SB() Roll(2.0*rollmain3) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
                   B(data,omega,l,phys_age);
       else if (nbuds == 4)
          produce  SB() Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB() 
                   SB() Roll(rollmain4) Down(1.0) Pitch(pitch)  B(data,omega+1,l1,phys_age) EB()
	           SB() Roll(2.0*rollmain4) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(3.0*rollmain4) Down(1.0) Pitch(pitch)  B(data,omega+1,l1,phys_age) EB()
                   B(data,omega,l,phys_age);
       else if (nbuds == 5)
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1,phys_age) EB() 
                   SB() Roll(rollmain5) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(2.0*rollmain5) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(3.0*rollmain5) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(4.0*rollmain5) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
                   B(data,omega,l,phys_age);
       else if (nbuds == 6)
          produce  SB() Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB() 
                   SB() Roll(rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
	           SB() Roll(2.0*rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(3.0*rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(4.0*rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(5.0*rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
                   B(data,omega,l,phys_age);
       else if (nbuds == 7)
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1,phys_age) EB() 
                   SB() Roll(rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(2.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(3.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(4.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(5.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(6.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
                   B(data,omega,l,phys_age);
       else if (nbuds >= 8)
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1,phys_age) EB() 
                   SB() Roll(rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(2.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(3.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(4.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(5.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
		   SB() Roll(6.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
                   SB() Roll(7.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1,phys_age) EB()
                   B(data,omega,l,phys_age);
    }//omega == 1
    else if (omega >= 2 && omega <= 100){
        //For branches start increasing physiological age
        data.phys_age += 1;
        phys_age += 1;
        if (nbuds == 0)
           produce B(data,omega,l,phys_age);
        else if ((is_architecture_change == true) && (phys_age >= architecture_change_year))
           //GROWTH ARCHITECTURE CHANGE IN branches: always only two new branches at the same plane to the opposite sides
           //NOTE: this can go forever unless growth conditions and functions set the segment length to 0 
           produce  SB() Roll(roll)   Turn(turn_branch)Down(3.0)   B(data,omega+1,l1,phys_age) EB() 
                    SB() Roll(2*roll)  Turn(turn_branch)Down(3.0)   B(data,omega+1,l1,phys_age) EB()
                    Down(-1) B(data,omega,l,phys_age);
        else if (nbuds == 1)
           produce  SB() Roll(roll)  Turn(turn_branch)Down(3.0)   B(data,omega+1,l1,phys_age) EB()
                    Down(-1) B(data,omega,l,phys_age);
        else if (nbuds == 2)
           //Two new branches at the same plane to the opposite sides
           produce  SB() Roll(roll)   Turn(turn_branch)Down(3.0)   B(data,omega+1,l1,phys_age) EB() 
                    SB() Roll(2*roll)  Turn(turn_branch)Down(3.0)   B(data,omega+1,l1,phys_age) EB()
                    Down(-1) B(data,omega,l,phys_age);
        else if (nbuds == 3)
           produce  SB() Roll(roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1,phys_age) EB() 
                    SB() Roll(2*roll)  Turn(turn_branch)Down(3.0)  B(data,omega+1,l1,phys_age) EB()
   		    SB() Roll(3.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1,phys_age) EB()
                    Down(-1) B(data,omega,l,phys_age);
        else if (nbuds == 4)
           produce  SB() Roll(roll)  Turn(turn_branch)Down(3.0)    B(data,omega+1,l1,phys_age) EB() 
                    SB() Roll(2.0*roll)  Turn(turn_branch)Down(3.0)    B(data,omega+1,l1,phys_age) EB()
   		    SB() Roll(3.0*roll)  Turn(turn_branch)Down(3.0)    B(data,omega+1,l1,phys_age) EB()
   		    SB() Roll(4.0*roll)  Turn(turn_branch)Down(3.0)     B(data,omega+1,l1,phys_age) EB()
                    Down(-1) B(data,omega,l,phys_age);
        else if (nbuds >= 5){
           produce  SB()  Roll(roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1,phys_age) EB() 
                    SB()  Roll(2.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1,phys_age) EB()
   		    SB()  Roll(3.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1,phys_age) EB()
   		    SB()  Roll(4.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1,phys_age) EB()
   		    SB()  Roll(5.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1,phys_age) EB()
                    Down(-1) B(data,omega,l,phys_age);
        }
    }//2 <= omega <= 10		    
    else{//omega>8
       //No more new branches for Gravelius order > 6
       data.phys_age += 1;
       phys_age += 1;
       produce B(data,omega,l,phys_age);
    }
 }//mode ==1
}//bud


//bend down axes branching from  the main  axis (omega == 2) until 90 degrees
//by producing a sequence of pitches 
Down(s):
{   
    if (s == -1){
      produce HDir(hdir) Down(1);
    }
    else if (s*pitch_down + pitch < max_pitch_down){
      produce Down(s+1.0)  Pitch(pitch_down);
    }
    else
      produce;
}


//This interpretation  is a ToDo for me.
//Now you must define at least one dummy interpretation.
//I must study how to alter the l2c to avoid  interpretation. 
//Interpretation separates model/symbols from 
//their (graphics) interpretation but for us LIGNUM is visualized by 
//Mika.  
interpretation:

F(s) :
{
  ;
}

//Closes the namespace (compiles simply to '}')
close Pine;
