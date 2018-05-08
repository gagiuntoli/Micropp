#include <stdlib.h> 
#include <iostream>
#include "micro.h"

using namespace std;

static Problem* micro = NULL;

extern "C"
{
  void micro_construct_(int* dim, int* size, int *cg_its, double *cg_tol)
  {
    micro = new Problem(*dim, size, *cg_its, *cg_tol); 
  }

  void micro_loc_hom_stress_(double *MacroStrain, double *MacroStress)
  {
    micro->loc_hom_Stress(MacroStrain, MacroStress); 
  }

  void micro_loc_hom_ctan_(double *MacroStrain, double *MacroCtan)
  {
    micro->loc_hom_Ctan(MacroStrain, MacroCtan); 
  }

  void micro_write_vtu_(int *time_step, int *elem)
  {
    micro->writeVtu(*time_step, *elem); 
  }
} 
