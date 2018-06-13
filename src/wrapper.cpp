#include <stdlib.h> 
#include <iostream>
#include "micro.h"

using namespace std;

static Problem* micro = NULL;

extern "C"
{
  void micro_construct_(int *dim, int size[3], int *micro_type, double *micro_params, int *mat_types, double *params)
  {
    micro = new Problem(*dim, size, *micro_type, micro_params, mat_types, params); 
  }

  void micro_loc_hom_stress_(int *macroGp_id, double *MacroStrain, double *MacroStress)
  {
    micro->loc_hom_Stress(*macroGp_id, MacroStrain, MacroStress); 
  }

  void micro_loc_hom_ctan_(int *macroGp_id, double *MacroStrain, double *MacroCtan)
  {
    micro->loc_hom_Ctan(*macroGp_id, MacroStrain, MacroCtan); 
  }

  void micro_output_(int *time_step, int *elem, int *macroGp_id, double *MacroStrain)
  {
    micro->output (*time_step, *elem, *macroGp_id, MacroStrain);
  }

  void micro_update_int_vars_(void)
  {
    micro->updateIntVars ();
  }

  void micro_get_non_linear_flag_(int *macroGp_id, int *non_linear)
  {
    micro->getNonLinearFlag (*macroGp_id, non_linear);
  }
} 
