#ifndef CROWNDENSITYGLOBALS_H
#define CROWNDENSITYGLOBALS_H
///\file LignumForestGlobals.h
///Global variables in in L-systems and in LignumForest.\footnote{This file gives examples of Doxygen grouping}
///
///To cut mutual dependency with CrownDensity global
///variables created in CrownDensity are moved to LignumForest
///in LignumForest namespace,
///\sa globalvariables.cc

#include <LGMSymbols.h>
#include <Point.h>
#include <ParametricCurve.h>
#include <Firmament.h>
#include <VoxelSpace.h>


///\brief Pine namespace
///
///L-system global variables needed in communication between L-system and Lignum tree,
///hence in the Pine namespace. See for example pine-em98.L and in its derivatives
///like pine-em98-branch-C.L
namespace Pine{
  ///\ingroup LIGNUMFORESTGLOBALGROUP
  ///\defgroup A_LSYSTEMGLOBALGROUP Global variables in L-system
  ///L-system global variables needed in communication between L-system and Lignum tree.
  ///These variables are in Pine namespace.
  ///@{
  ///\brief Tree age
  extern double L_age;
  ///\brief Tree height
  extern double L_H;
  ///\brief A guess for segment shortening
  extern const double l1;
  ///\brief Number of buds as a function of foliage mass
  extern ParametricCurve fnbuds;
  ///\brief Adjust number of buds
  ///Adjust the number of buds returned by \p fnbuds with relative light
  extern ParametricCurve fnbudslight;
  ///\brief Number of buds as function of relative light
  ///\note Used the reset number of buds to zero in L system if f(ip) = 0.
  extern ParametricCurve fipbud;
  ///\brief Set architecture change on/off
  extern bool is_architecture_change;
  ///\brief Set architecture change year
  extern int architecture_change_year;
  ///@}
}


///\brief LignumForest namespace
///
///Global variables used in LignumForest initialized by GrowthLoop.
///To cut mutual dependency with CrownDensity global
///variables from CrownDensity were moved to LignumForest
///and explicitely referred with LignumForest namespace.
namespace LignumForest{
  ///\defgroup LIGNUMFORESTGLOBALGROUP Global variables
  ///Global variables for various experiments to model and adjust forest growth
  ///@{
  ///@}
  ///\ingroup LIGNUMFORESTGLOBALGROUP 
  ///\defgroup B_HEIGHTVARS Global variables for initial sapling height
  ///@{
  ///\brief Initial height ot trees
  extern double H_0_ini;
  ///\brief Variation of initial tree height
  extern double H_var_ini;
  ///@}
  ///\ingroup LIGNUMFORESTGLOBALGROUP 
  ///\defgroup C_INITIALSAPLINGBUDVARS Variables for buds in the initial saplings
  ///@{
  ///\brief Set bud variation on or off \sa rel_bud
  extern bool bud_variation;
  ///\brief The effect of crowding on buds
  ///
  ///The effect of crowding on number of lateral
  ///buds via function bud_view_f.
  ///\pre Requires LignumForest::bud_variation
  ///\sa LignumForest::bud_variation
  extern double rel_bud;
  ///\brief Variation in the initial max number of buds.
  extern int n_buds_ini_max;
  ///\brief Variation in the initial min number of buds.
  extern int n_buds_ini_min;
  ///\deprecated Variation in branch angle. Set in \sa pine-em98.L
  extern double branch_angle;
  ///@}
  
}

namespace LignumForest{
  ///\ingroup LIGNUMFORESTGLOBALGROUP
  ///\defgroup D_RANDOMVARS Variables for randomness
  ///@{
  ///\brief Seed for ran3 function
  extern int ran3_seed;
  ///\brief Random component in segment length.
  ///
  ///Use random component to apply conditionally any function implementing segment elongation
  extern bool is_random_length;
  ///\brief Ad hoc lengthening experiment of shoots at crown base.
  ///\sa CrownDensity::adoc
  extern bool is_adhoc;
  ///\brief Function for ad hoc lengthening experiment of shoots at crown base.
  ///\sa LignumForest::is_adhoc
  extern ParametricCurve adhoc;
  ///@}
  
  ///\brief Global variable for EBH model is needed e.g. in bybranches.cc-
  ///
  ///\todo Not (yet) used in LignumForest
  ///\note GrowthLoop.h has growthloop_ebh_mode as a class data member
  ///\note GrowthLoopI.h sets `ebh_mode` in when parsing command line
  extern int ebh_mode;       
  ///\brief Radiation use efficiency experiment
  ///\attention Not (yet) used in LignumForest
  extern LGMdouble max_rueqin;
  ///\ingroup LIGNUMFORESTGLOBALGROUP
  ///\defgroup E_RELATIVELIGHT Segment length functions
  ///@{
  ///\brief Relative light and Gravelius order functions
  ///
  ///\attention The Fip, effect of relative light, and effect of Gravielius order
  ///on segment length are in use but are read from MetaFiles.
  ///\brief Fip function, effect of relative light on segment
  ///\deprecated Fip functions are read from MetaFiles 
  extern ParametricCurve fip_mode;
  ///\brief Gravelius order function f(go) effecting segment length after mode change
  ///\deprecated Gravelius order functions are read from MetaFiles
  extern ParametricCurve fgo_mode;
  ///@}
  ///\ingroup LIGNUMFORESTGLOBALGROUP
  ///\defgroup F_GROWTHMODE Growth mode change variables
  ///@{
  ///\brief Reinitialize trees
  ///
  ///Reread MetaFiles and reinitialize trees.
  ///\brief Set mode change on/off
  ///\sa GrowthLoop::insertModeChangeYears()
  extern bool is_mode_change;
  ///@}
  ///\ingroup LIGNUMFORESTGLOBALGROUP
  ///\defgroup G_GROWTHSPACE Growth space experiment variables
  ///@{
  ///\brief Additional Firmament for  VoxelSpace `LignumForest::space_occupancy`.
  ///\sa LignumForest::space_occupancy
  extern Firmament dummy_firm;
  ///\brief Additional VoxelSpace for space occupancy
  ///\sa sky::dummy_firm
  extern VoxelSpace space_occupancy;
  ///\brief Three options \p space0, \p space1 and \p space2, for three space occupancy model experiments
  ///\name Space occupancy models
  ///@{
  extern bool space0;
  extern bool space1;
  extern bool space2;
  ///\brief Look ahead distance for a bud  in \p space2 model
  extern double space2_distance;
  ///@}
  ///@}
  ///\ingroup  LIGNUMFORESTGLOBALGROUP
  ///\defgroup H_STEMSHAPE Crown and stem shape experiment
  ///@{
  ///\brief Height crown base
  extern double global_hcb;
  ///\brief For change in base diameter, for the command line option -heightFun
  ///\attention The class GrowthLoop has dDb as data member 
  extern double dDb;
  ///@}

  ///\ingroup LIGNUMFORESTGLOBALGROUP
  ///\defgroup Z_OBSELETE Variables not in use
  ///@{
  ///\deprecated Not used in LignumForest
  ///\todo Renmove from LignumForest
  extern bool is_height_function;
  ///@}
  
  ///\ingroup LIGNUMFORESTGLOBALGROUP
  ///\defgroup I_BUDVIEW Bud lookahead view experiment
  ///@{
  ///\brief Bud view function
  ///
  ///Global variable to convey the Bud View Function to L-system
  ///and if LignumForest::bud_view_f is in use
  ///\attention initialized with "bvf.fun" in globalvariables.cc
  extern ParametricCurve bud_view_f;
  ///\brief Boolean flag to set bud view function on or off.
  ///\sa bud_vew_f
  extern bool is_bud_view_function;
  ///\brief The angle of branching after architecture change. Set in L system.
  extern double max_turn_in_architecture_change;
  ///@}
  ///\ingroup LIGNUMFORESTGLOBALGROUP
  ///\defgroup J_BUTTSWELL Butt swell experiment
  ///@{
  ///\brief Butt swell model in diameter growth.
  ///
  ///Adjust \f$LGPq\f$ in butt swell model with a coefficient \f$c\f$: \f$LGPq = c*LGPq\f$ and \f$0 \le c \le 1\f$.
  ///The default value for \f$c\f$ is 1.0 (i.e. no effect).
  ///The butt swell model is implemented in ScotsPineSegment. The most straigthforward way
  ///to implement adjustment parameters is to use global variables in LignumForest namespace.
  ///\note Butt swell is in use
  ///\sa ScotsPineSegment::aging() GrowthLoop::usage() GrowthLoop::parseCommanLine() LignumForest::butt_swell_start 
  extern double butt_swell_coeff;
  ///\brief Tree age to start butt swell, default value is INT_MAX (i.e. never).
  ///\sa ScotsPineSegment::aging() GrowthLoop::usage()  GrowthLoop::parseCommanLine() LignumForest::butt_swell_coeff
  extern int butt_swell_start;
  ///@}
  ///\ingroup LIGNUMFORESTGLOBALGROUP
  ///\defgroup K_BUDTERMINATION 11. Terminate buds
  ///Terminate buds grown out of VoxelSpace
  ///@{
  ///\brief Terminate buds
  ///
  ///Boolean flag. Terminate buds grown out of VoxelSpace, default False.
  ///The command line option `-terminate_buds` sets the value to True.
  ///\note Bud termination is in use
  extern bool terminate_escaped_buds;
  ///@}

  ///Limiting length of new shoots
  extern double g1maxL, g2maxL;
  extern int length_limit_year;

}//end namespace
#endif
