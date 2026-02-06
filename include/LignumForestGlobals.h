#ifndef CROWNDENSITYGLOBALS_H
#define CROWNDENSITYGLOBALS_H
/// \file LignumForestGlobals.h
/// \brief Global variables in in L-systems and in LignumForest.
///
///\note To remove mutual dependency with CrownDensity global
///variables created in CrownDensity are moved to LignumForest
///in LignumForest namespace,
///\sa globalvariables.cc

#include <LGMSymbols.h>
#include <Point.h>
#include <ParametricCurve.h>
#include <Firmament.h>
#include <VoxelSpace.h>
using namespace voxelspace;

//This seems to be a convinient way to group items.
//First, define groups and subgroups with brief descriptions.
//Then add items to groups with *addtogroup* (multiple items)
//or *ingroup* (single item).
//
//The groups will appear in the order they appear for Doxygen
//or Doxygen can sort the groups (SORT_GROUP_NAMES = YES).
/// \defgroup GLOBALGROUP Global variables
/// \brief Global variables in LignumForest.
///
///Global variables passing information between LignumForest and L-systems
///and for various experiments to model and adjust forest growth.
/// @{
/// \defgroup GROWTHMODE Growth mode
/// \brief Variables for growth mode change experiment
/// \defgroup GROWTHSPACE Growth space
/// \brief Variables for growth space experiment
///
/// Three growth space occupancy models.
///
/// \defgroup BUDSPACE Bud space
/// \brief Variables for the bud growth sector clearance experiment
/// \defgroup BUDTERMINATION Bud termination
/// \brief Variables to cease buds not in VoxelSpace
///
/// Buds grown outside the \e original VoxelSpace are marked Lignum::DEAD.
///
/// \defgroup BUTTSWELL Butt swell
/// \brief Variables for butt swell model
/// \defgroup SHAPE Crown and stem shape
/// \brief Variables for crown and stem shape experiments
/// \defgroup INITIALSAPLINGS Initial saplings
/// \brief Variables for initial saplings
///
/// Initial sapling (i.e. tree) length and length variation.
///
/// \defgroup INITIALBUDS Initial buds
/// \brief Variables for buds in the initials saplings
///
/// Variables for minimum and maximum number of buds, bud effect on crowding, brancing angle
/// and variation in branching angle.
///
/// \defgroup LSYSTEM  L-system
/// \brief Variables for L-system
///
/// L-system global variables needed in communication between L-system and Lignum tree.
/// These variables are in Pine namespace. See for example \e pine-em98-branch-C.L.
///
/// \defgroup RANDOM Randomness
/// \brief Variables for randomness
///
/// \defgroup ALLOCATION Resource allocation
/// \brief Variables for resource allocation experiments
///
/// Extended Borchert–Honda and Radiation use efficiency models.
///
/// \defgroup SEGMENTLIMIT Segment length limit
/// \brief Variables limiting growth in new segments experiment
///
/// Impose max length for first and second order branch segments.
/// First branch segments and second branch segments can have different max lengths.
/// Main axis and higher branches are \e not affected. The model behaviour is defined
/// in the file LignumForest::SEGMENT_LENGTH_LIMIT_FILE.
/// \attention Be careful not to set too low values of max segment lengths.
/// Undesired effects may occur through low foliage mass.
///
/// \defgroup DEPRECATED Deprecated variables
/// \brief Variables not used
/// \attention The Fip, effect of relative light, and effect of Gravielius order on segment length
/// are read from MetaFiles.
/// @}

///\brief L-system namespace
namespace Pine{
  /// \addtogroup LSYSTEM
  /// @{
  ///\brief Tree age
  extern double L_age;
  ///\brief Tree height
  extern double L_H;
  ///\brief A guess for segment shortening
  extern const double l1;
  ///\brief Number of buds as a function of foliage mass
  extern ParametricCurve fnbuds;
  ///\brief Adjust number of buds
  ///
  ///Adjust the number of buds returned by \p LignumForest::fnbuds with relative light
  extern ParametricCurve fnbudslight;
  ///\brief Number of buds as function of relative light
  ///
  ///\note Used the reset number of buds to zero in L system if f(ip) = 0.
  extern ParametricCurve fipbud;
  ///\brief Set architecture change on/off
  extern bool is_architecture_change;
  ///\brief Set architecture change year
  extern int architecture_change_year;
  /// @}
}


///\brief LignumForest namespace
///
///Global variables used in LignumForest initialized by GrowthLoop.
///
///To cut mutual dependency with CrownDensity global
///variables from CrownDensity are moved to LignumForest
///and explicitely referred with LignumForest namespace.
namespace LignumForest{
  ///\addtogroup INITIALSAPLINGS 
  ///@{
  ///\brief Initial height ot trees
  extern double H_0_ini;
  ///\brief Variation of initial tree height
  extern double H_var_ini;
  ///@}
  ///\addtogroup INITIALBUDS
  ///@{
  ///\brief Set bud variation on or off \sa LignumForest::rel_bud
  extern bool bud_variation;
  ///\brief The effect of crowding on buds
  ///
  ///The effect of crowding on number of lateral
  ///buds via function LignumForest::bud_view_f.
  ///\pre Requires LignumForest::bud_variation
  ///\sa LignumForest::bud_variation
  extern double rel_bud;
  ///\brief Variation in the initial max number of buds.
  extern int n_buds_ini_max;
  ///\brief Variation in the initial min number of buds.
  extern int n_buds_ini_min;
  ///\deprecated Variation in branch angle. Set in \sa pine-em98.L
  extern double branch_angle;
  ///\brief The angle of branching after architecture change. Used in L system.
  extern double max_turn_in_architecture_change;
  ///@}
  
}

namespace LignumForest{
  ///\addtogroup RANDOM
  ///@{
  ///\brief Seed for `ran3` function
  extern int ran3_seed;
  ///\brief Random component in segment length.
  ///
  ///Use random component to apply conditionally any function implementing segment elongation
  extern bool is_random_length;
  ///\brief Ad hoc lengthening experiment of shoots at crown base.
  ///\sa CrownDensity::adoc
  extern bool is_adhoc;
  ///\brief Function for *ad hoc* lengthening experiment of shoots at crown base.
  ///\sa LignumForest::is_adhoc
  extern ParametricCurve adhoc;
  ///@}
  ///\addtogroup ALLOCATION
  ///@{
  ///
  ///\brief EBH model is used e.g. in bybranches.cc
  ///\todo Not (yet) used in LignumForest
  ///\note GrowthLoop.h has growthloop_ebh_mode as a class data member
  ///\note GrowthLoopI.h sets `ebh_mode` in when parsing command line
  extern int ebh_mode;       
  ///\brief Radiation use efficiency experiment
  ///\attention Not (yet) used in LignumForest
  extern LGMdouble max_rueqin;
  ///@}
  ///\addtogroup GROWTHMODE
  ///@{
  ///\brief Reinitialize trees
  ///
  ///Set mode change on/off to rereade
  ///MetaFiles and reinitialize trees.
  ///\sa GrowthLoop::insertModeChangeYears()
  extern bool is_mode_change;
  ///@}
  ///\addtogroup GROWTHSPACE
  ///@{
  ///\brief Additional Firmament for  VoxelSpace `LignumForest::space_occupancy`.
  ///\sa LignumForest::space_occupancy
  extern Firmament dummy_firm;
  ///\brief Additional VoxelSpace for space occupancy
  ///\sa sky::dummy_firm
  extern VoxelSpace space_occupancy;
  ///\name Space occupancy models
  ///\brief Select one of the three options \p space0, \p space1 or \p space2 for the space occupancy model experiment.
  ///@{
  extern bool space0;
  extern bool space1;
  extern bool space2;
  ///\brief Look ahead distance for a bud  in \p LignumForest::space2 model
  extern double space2_distance;
  ///@}
  ///@}
  ///\addtogroup SHAPE
  ///@{
  ///\brief Height crown base
  extern double global_hcb;
  ///\brief For change in base diameter, for the command line option -heightFun
  ///\attention The class GrowthLoop has dDb as data member 
  extern double dDb;
  ///@}
  ///\addtogroup DEPRECATED 
  ///@{
  ///\brief Height functionn experiment in use or not
  ///\deprecated Not used in LignumForest
  ///\todo Remove from LignumForest
  extern bool is_height_function;
  ///@}
  ///\addtogroup DEPRECATED 
  ///@{
  ///\brief Fip function
  ///
  ///Effect of relative light on segment length
  ///\deprecated Fip functions are read from MetaFiles 
  extern ParametricCurve fip_mode;
  ///\brief Gravelius order function
  ///
  ///Gravelius order effect on segment length 
  ///\deprecated Gravelius order functions are read from MetaFiles
  extern ParametricCurve fgo_mode;
  ///@}
  ///\addtogroup BUDSPACE
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
  ///@}
  ///\addtogroup BUTTSWELL
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
  ///\addtogroup BUDTERMINATION
  ///@{
  ///\brief Variable for command line 
  ///
  ///Boolean flag. Terminate buds grown out of the *original* VoxelSpace, default False.
  ///The command line option `-terminate_buds` sets the value to *true*.
  ///\sa GrowthLoop<TREE,TS,BUD,LSYSTEM>::terminateEscapedBuds().
  extern bool terminate_escaped_buds;
  ///@}
  ///\addtogroup SEGMENTLIMIT 
  ///@{
  ///\brief File defining max segment length model behaviour.
  ///
  ///File sets values for LignumForest::g1maxL, LignumForest::g2maxL
  ///and LignumForest::length_limit_year.
  const string SEGMENT_LENGTH_LIMIT_FILE="dhlimit.txt";
  ///First branch segments max length. Use unrealistic high value to have no effect.
  extern double g1maxL;
  ///Second branch segments max length. Use unrealistic high value to have no effect.
  extern double g2maxL;
  ///\brief Max segment length start year.
  ///
  ///Start year the max segment lengths LignumForest::g1maxL and LignumForest::g2maxL take effect.
  ///Value larger than the simulation time annuls the max segment length model calculations.
  extern int length_limit_year;
  ///@}
}//end namespace
#endif
