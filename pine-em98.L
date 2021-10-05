#include <iostream>
using namespace std;
//Defines ran3 etc
#include <mathsym.h>
using namespace cxxadt;
#include <LGMUnits.h>
using namespace Lignum;
#include <Pine.h>
//Defines headers necessary for lsystem
#include <lengine.h>

//These globals transfer information between L-system
//and other part of program. 
double L_age, L_H;

extern int ran3_seed;

extern double H_0_ini, H_var_ini;          //For variation of initial heights and
extern int n_buds_ini_min, n_buds_ini_max;   // and number of buds (in .L file)
extern double  rel_bud;
extern bool  bud_variation;
extern double branch_angle;

//If foliage area density in the view cone of the bud affects
//the number of lateral buds created
extern ParametricCurve bud_view_f;
extern bool is_bud_view_function;   ///if it is in use 



//============

//Start the namespace (compiles to "namespace Pine{")
//The matching "close" statement closes the namespace 
open Pine;

//Defines the lsystem rewriting process. 
#include <lsys.h>

//a guess for the segment shortening
const double l1 = 0.9;


const double PI = 3.1415926535897932384;
//const double pitch = 45.0*PI/180.0; //the angle of branching in the main axis    = 45
double pitch = 45.0*PI/180.0; //the angle of branching in the main axis    = 45
//double pitch = branch_angle;    //branch_agle may vary, it is set in GrowthLoopI.h
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
//NOTE for pine-forest: As the 'mode' is a global variable, it is visible to all L-systems.
//Consequently each l[i].derive() will change the mode and alternately the L-systems 
//will receive 0 and 1 as the mode. This is not what we want but what is needed is an additional
//explicit call to l[i].endEach() after each l[i].derive() in pine-forest.cc that resets the mode.
int mode = 0; //0 = new segment, 1 = new buds

//Number of buds as the function of foliage mass
ParametricCurve fnbuds("fnbuds.fun");

//Adjust the number of buds, more light more buds, less light less buds
ParametricCurve fnbudslight("fnbudslight.fun");
//Fip for buds, the death zone
ParametricCurve fipbud("fip-bud.fun");
module F(double);
//Bud(data, omega, initial length)
module B(PineBudData,int,double);
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

pitch = branch_angle;

  PineBudData data(ALIVE,0.0,1.0,0.1);
  mode = 0;
  int n_buds_ini = n_buds_ini_min + (int)(ran3(&ran3_seed)*(double)(n_buds_ini_max-n_buds_ini_min+1));

  if(n_buds_ini <= 1)
        produce F(H_0_ini+ran3(&ran3_seed)*H_var_ini) SB() Down(1.0) Pitch(pitch) B(data,2,l1) EB() 
                  B(data,1,1.0); 
  else if(n_buds_ini == 2)
        produce F(H_0_ini+ran3(&ran3_seed)*H_var_ini) SB() Down(1.0) Pitch(pitch) B(data,2,l1) EB() 
                  SB() Roll(rollmain2) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
                  B(data,1,1.0); 

  else if(n_buds_ini == 3)
       produce F(H_0_ini+ran3(&ran3_seed)*H_var_ini) SB() Down(1.0) Pitch(pitch) B(data,2,l1) EB() 
                  SB() Roll(rollmain3) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
	          SB() Roll(2.0*rollmain3) Down(1.0) Pitch(pitch) B(data,2,1.0) EB()
                  B(data,1,1.0); 

  else if(n_buds_ini == 4)
       produce F(H_0_ini+ran3(&ran3_seed)*H_var_ini) SB() Down(1.0) Pitch(pitch) B(data,2,l1) EB() 
                  SB() Roll(rollmain4) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
	          SB() Roll(2.0*rollmain4) Down(1.0) Pitch(pitch) B(data,2,1.0) EB()
		  SB() Roll(3.0*rollmain4) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
                  B(data,1,1.0); 
  else if(n_buds_ini == 5)
       produce F(H_0_ini+ran3(&ran3_seed)*H_var_ini) SB() Down(1.0) Pitch(pitch) B(data,2,l1) EB() 
                  SB() Roll(rollmain5) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
	          SB() Roll(2.0*rollmain5) Down(1.0) Pitch(pitch) B(data,2,1.0) EB()
		  SB() Roll(3.0*rollmain5) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
		  SB() Roll(4.0*rollmain5) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
                  B(data,1,1.0); 
  else if(n_buds_ini == 6)
       produce F(H_0_ini+ran3(&ran3_seed)*H_var_ini) SB() Down(1.0) Pitch(pitch) B(data,2,l1) EB() 
                  SB() Roll(rollmain6) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
	          SB() Roll(2.0*rollmain6) Down(1.0) Pitch(pitch) B(data,2,1.0) EB()
		  SB() Roll(3.0*rollmain6) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
		  SB() Roll(4.0*rollmain6) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
		  SB() Roll(5.0*rollmain6) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
                  B(data,1,1.0); 

  else if(n_buds_ini == 7)
       produce F(H_0_ini+ran3(&ran3_seed)*H_var_ini) SB() Down(1.0) Pitch(pitch) B(data,2,l1) EB() 
                  SB() Roll(rollmain7) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
	          SB() Roll(2.0*rollmain7) Down(1.0) Pitch(pitch) B(data,2,1.0) EB()
		  SB() Roll(3.0*rollmain7) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
		  SB() Roll(4.0*rollmain7) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
 		  SB() Roll(5.0*rollmain7) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
 		  SB() Roll(6.0*rollmain7) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB() 
                 B(data,1,1.0); 

  else
  produce F(H_0_ini+ran3(&ran3_seed)*H_var_ini) SB() Down(1.0) Pitch(pitch) B(data,2,l1) EB() 
                  SB() Roll(rollmain8) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
	          SB() Roll(2.0*rollmain8) Down(1.0) Pitch(pitch) B(data,2,1.0) EB()
		  SB() Roll(3.0*rollmain8) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
		  SB() Roll(4.0*rollmain8) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
		  SB() Roll(5.0*rollmain8) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
		  SB() Roll(6.0*rollmain8) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
		  SB() Roll(7.0*rollmain8) Down(1.0) Pitch(pitch)  B(data,2,1.0) EB()
                  B(data,1,1.0);
}

StartEach:
{
}

EndEach:
{
//mode = (mode+1) % 2;
} 

B(data,omega,l):
{
pitch = branch_angle;
  if (mode == 0){
   if (data.state == ALIVE){
     //Note there are no start and end branches
     //These are genrerated in mode 1
     if (omega==1)
       produce F(l) Roll(roll) B(data,omega,l);
     else
       produce F(l) B(data,omega,l);
   }
   else{
       produce B(data,omega,l);
   }
  }
else{//mode == 1
//Here also the effect of crowding (in terms of foliage area density)
//on number of lateral buds via function bud_view_f

    double rel = 1.0;
    if(bud_variation) {
      rel += rel_bud*(ran3(&ran3_seed)-0.5)/0.5;
      if(rel < 0.0) rel = 0.0;
    }

    int nbuds = static_cast<int>(rel*fnbuds(data.fm));

    //then adjust the number of buds with relative light  
    //nbuds = static_cast<int>(max(nbuds+fnbudslight(data.ip),0.0));
    //cout << "ip " << data.ip << " fip(ip) " << fip(data.ip) 
   //     << " FM " << data.fm << " fnbuds.fun " << nbuds <<endl;

    if(fipbud(data.ip) == 0.0) nbuds = 0;

    if (nbuds == 1 || nbuds == 2){
       if (fipbud(data.ip) < 0.6){
           nbuds = 0;
       }
    }

    // Local density (foliage area density) possibly affects number of buds  
    if (is_bud_view_function) {
      nbuds = static_cast<int>(static_cast<double>(nbuds) * bud_view_f(data.view));
    }

    if (nbuds == 0 && L_age < 4)  
       nbuds = 1;                

    if (data.state == DEAD){
       //cout << "Already dead " << endl;
       produce B(data,omega,l);
    }
    else if (nbuds == 0){
       data.state = DEAD;
       //cout << "No new buds, DEAD" <<endl;
       produce B(data,omega,l);
    }
    else if (omega == 1){
       if (nbuds == 1)
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1) EB()
                   B(data,omega,l);
       else if (nbuds == 2)
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1) EB() 
                   SB() Roll(rollmain2) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
                   B(data,omega,l);
       else if (nbuds == 3)
          produce  SB() Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB() 
                   SB() Roll(rollmain3) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
	           SB() Roll(2.0*rollmain3) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
                   B(data,omega,l);
       else if (nbuds == 4)
          produce  SB() Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB() 
                   SB() Roll(rollmain4) Down(1.0) Pitch(pitch)  B(data,omega+1,l1) EB()
	           SB() Roll(2.0*rollmain4) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(3.0*rollmain4) Down(1.0) Pitch(pitch)  B(data,omega+1,l1) EB()
                   B(data,omega,l);
       else if (nbuds == 5){
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1) EB() 
                   SB() Roll(rollmain5) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(2.0*rollmain5) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(3.0*rollmain5) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(4.0*rollmain5) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
                   B(data,omega,l);
       }
       else if (nbuds == 6)
          produce  SB() Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB() 
                   SB() Roll(rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
	           SB() Roll(2.0*rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(3.0*rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(4.0*rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(5.0*rollmain6) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
                   B(data,omega,l);
       else if (nbuds == 7)
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1) EB() 
                   SB() Roll(rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(2.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(3.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(4.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(5.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(6.0*rollmain7) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
                   B(data,omega,l);
       else if (nbuds >= 8)
          produce  SB() Down(1.0) Pitch(pitch)B(data,omega+1,l1) EB() 
                   SB() Roll(rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(2.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(3.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(4.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(5.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
		   SB() Roll(6.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
                   SB() Roll(7.0*rollmain8) Down(1.0) Pitch(pitch) B(data,omega+1,l1) EB()
                   B(data,omega,l);
    }//omega == 1
    else if (omega >= 2 && omega <= 6){
        if (nbuds == 0)
           produce B(data,omega,l);
        else if (nbuds == 1)
           produce  SB() Roll(roll)  Turn(turn_branch)Down(3.0)   B(data,omega+1,l1) EB()
                    Down(-1) B(data,omega,l);
        else if (nbuds == 2)
           //two new branches at the same plane to the opposite sides
           produce  SB() Roll(roll)   Turn(turn_branch)Down(3.0)   B(data,omega+1,l1) EB() 
                    SB() Roll(2*roll)  Turn(turn_branch)Down(3.0)   B(data,omega+1,l1) EB()
                    Down(-1) B(data,omega,l);
        else if (nbuds == 3)
           produce  SB() Roll(roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB() 
                    SB() Roll(2*roll)  Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(3.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
                    Down(-1) B(data,omega,l);
        else if (nbuds == 4)
           produce  SB() Roll(roll)  Turn(turn_branch)Down(3.0)    B(data,omega+1,l1) EB() 
                    SB() Roll(2.0*roll)  Turn(turn_branch)Down(3.0)    B(data,omega+1,l1) EB()
   		    SB() Roll(3.0*roll)  Turn(turn_branch)Down(3.0)    B(data,omega+1,l1) EB()
   		    SB() Roll(4.0*roll)  Turn(turn_branch)Down(3.0)     B(data,omega+1,l1) EB()
                    Down(-1) B(data,omega,l);
        else if (nbuds == 5){
           produce  SB()  Roll(roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB() 
                    SB()  Roll(2.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB()  Roll(3.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB()  Roll(4.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB()  Roll(5.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
                    Down(-1) B(data,omega,l);
        }
        else if (nbuds == 6)
           produce  SB() Roll(roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB() 
                    SB() Roll(2.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(3.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(4.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(5.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
                    SB() Roll(6.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
                    Down(-1) B(data,omega,l);
        else if (nbuds == 7)
           produce  SB() Roll(roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB() 
                    SB() Roll(2.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(3.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(4.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(5.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
                    SB() Roll(6.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
		    SB() Roll(7.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
                    Down(-1) B(data,omega,l);
        else if (nbuds >= 8)
           produce  SB() Roll(roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB() 
                    SB() Roll(2.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(3.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(4.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
   		    SB() Roll(5.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
                    SB() Roll(6.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
		    SB() Roll(7.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
                    SB() Roll(8.0*roll)   Turn(turn_branch)Down(3.0)  B(data,omega+1,l1) EB()
                    Down(-1) B(data,omega,l);
    }//2 <= omega <= 6  		    
    else{
       produce B(data,omega,l);
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
