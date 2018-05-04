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

    if (micro->flag_print_wrapper == true) {
      cout<< "dim = " << *dim << endl; 
      cout<< "size = " << size[0] << " " << size[1] << " " << size[2] << endl; 
      cout<< "cg_its = " << *cg_its << endl; 
      cout<< "cg_tol = " << *cg_tol << endl; 
    }
  }

  void micro_loc_hom_stress_(double *MacroStrain, double *MacroStress)
  {
    micro->loc_hom_Stress(MacroStrain, MacroStress); 
  }

  void micro_loc_hom_ctan_(double *MacroStrain, double *MacroCtan)
  {
    micro->loc_hom_Ctan(MacroStrain, MacroCtan); 
  }
} 
