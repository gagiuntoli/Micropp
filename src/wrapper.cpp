/*
 *  MicroPP : finite element library to solve microstructure problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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

  void micro_loc_hom_ctan_linear_(double *Ctan)
  {
    micro->loc_hom_Ctan_Linear (Ctan);
  }

  void micro_loc_hom_stress_linear_(double *strain, double *stress)
  {
    micro->loc_hom_Stress_Linear (strain, stress);
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

  void micro_get_params_nr_(int *NR_its, double *NR_norm, int *NR_non_linear)
  {
    micro->getParams_NR (NR_its, NR_norm, NR_non_linear);
  }

  void micro_get_int_vars_(int *macroGp_id, int *n, int *int_vars)
  {
    micro->getIntVars (*macroGp_id, *n, int_vars);
  }

  void micro_get_max_ftrial_(double *ftrial_max)
  {
    micro->getMaxFtrial (ftrial_max);
  }

  void micro_get_max_invariant_(double *invariant_max)
  {
    micro->getMaxInvariant (invariant_max);
  }

  void micro_update_ctan_static_(void)
  {
    micro->updateCtanStatic ();
  }

  void micro_get_ctan_static_(int *MacroGp_id, double *Ctan)
  {
    micro->getCtanStatic (*MacroGp_id, Ctan);
  }

} 
