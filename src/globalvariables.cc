#include <mathsym.h>
#include <LignumForestGlobals.h>
///\file globalvariables.cc

namespace Pine{
  double L_age = 0.0;
  double L_H = 0.0;
  const double l1 = 0.9;
  //ParametricCurve fnbuds;
  //ParametricCurve fnbudslight;
  //ParametricCurve fipbud;
  bool is_architecture_change = false;
  int architecture_change_year = INT_MAX;
}


namespace LignumForest{
  double H_0_ini = 0.0;
  double H_var_ini = 0.0;
  int n_buds_ini_max = 0.0;
  int n_buds_ini_min = 0.0; 
  bool bud_variation = false;
  double rel_bud = 0.0;
  double branch_angle = 45.0 * PI_VALUE / 180.0;
  double max_turn_in_architecture_change = 80.0*PI_VALUE/180.0;
  double butt_swell_coeff = 1.0;
  int butt_swell_start = INT_MAX;
  double g1maxL = 1.0, g2maxL = 1.0;
  int length_limit_year = 1000;      

}

namespace LignumForest{
  int ran3_seed=-1234567;
  bool is_random_length = false;
  bool is_adhoc = false;
  ParametricCurve adhoc("adhoc.fun");
  int ebh_mode = 0;          
  LGMdouble max_rueqin;
  bool is_mode_change = false;
  Firmament dummy_firmament;
  VoxelSpace space_occupancy(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
			     0.1,0.1,0.1,5,5,5,dummy_firmament);
  double global_hcb;
  double dDb;
  bool space0 = false;
  bool space1 = false;
  bool space2 = false;
  double space2_distance = 0.3;
  bool is_height_function = false;
  ParametricCurve bud_view_f("bvf.fun");
  bool is_bud_view_function = false;
  bool terminate_escaped_buds = false;
}//end namespace
