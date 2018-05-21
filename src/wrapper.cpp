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

  void micro_loc_hom_stress_(double *MacroStrain, double *MacroStress)
  {
    micro->loc_hom_Stress(MacroStrain, MacroStress); 
  }

  void micro_loc_hom_ctan_(double *MacroStrain, double *MacroCtan)
  {
    micro->loc_hom_Ctan(MacroStrain, MacroCtan); 
  }

  void micro_output_(int *time_step, int *elem, int *macro_gp_global, double *MacroStrain)
  {
    micro->output (*time_step, *elem, *macro_gp_global, MacroStrain);
  }
} 
