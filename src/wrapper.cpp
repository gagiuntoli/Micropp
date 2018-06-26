/*
 *  MicroPP : 
 *  Finite element library to solve microstructural problems for composite materials.
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
  void micropp_construct_(int *dim, int size[3], int *micro_type, double *micro_params, int *mat_types, double *params)
  {
    micro = new Problem(*dim, size, *micro_type, micro_params, mat_types, params); 
  }

  void micropp_loc_hom_stress_(int *Gauss_ID, double *MacroStrain, double *MacroStress)
  {
    micro->loc_hom_Stress(*Gauss_ID, MacroStrain, MacroStress); 
  }

  void micropp_loc_hom_ctan_(int *Gauss_ID, double *MacroStrain, double *MacroCtan)
  {
    micro->loc_hom_Ctan(*Gauss_ID, MacroStrain, MacroCtan); 
  }

  void micropp_loc_hom_ctan_linear_(double *Ctan)
  {
    micro->loc_hom_Ctan_Linear (Ctan);
  }

  void micropp_loc_hom_stress_linear_(double *strain, double *stress)
  {
    micro->loc_hom_Stress_Linear (strain, stress);
  }

  void micropp_output_(int *time_step, int *Gauss_ID, double *MacroStrain)
  {
    micro->output (*time_step, *Gauss_ID, MacroStrain);
  }

  void micropp_get_non_linear_flag_(int *Gauss_ID, int *non_linear)
  {
    micro->getNonLinearFlag (*Gauss_ID, non_linear);
  }

  void micropp_get_params_nr_(int *NR_its, double *NR_norm, int *NR_non_linear)
  {
    micro->getParams_NR (NR_its, NR_norm, NR_non_linear);
  }

  void micropp_get_int_vars_(int *Gauss_ID, int *n, int *int_vars)
  {
    micro->getIntVars (*Gauss_ID, *n, int_vars);
  }

  void micropp_update_ctan_static_(void)
  {
    micro->updateCtanStatic ();
  }

  void micropp_get_ctan_static_(int *Gauss_ID, double *Ctan)
  {
    micro->getCtanStatic (*Gauss_ID, Ctan);
  }

  void micropp_set_macro_strain_(int *Gauss_ID, double *MacroStrain)
  {
    micro->setMacroStrain (*Gauss_ID, MacroStrain);
  }

  void micropp_localize_homogenize_(void)
  {
    micro->localizeHomogenize();
  }

  void micropp_get_macro_stress_(int *Gauss_ID, double *MacroStress)
  {
    micro->getMacroStress(*Gauss_ID, MacroStress);
  }

  void micropp_get_macro_ctan_(int *Gauss_ID, double *MacroCtan)
  {
    micro->getMacroCtan(*Gauss_ID, MacroCtan);
  }

  void micropp_update_internal_variables_(void)
  {
    micro->updateInternalVariables ();
  }

  void micropp_update_int_vars_(void)
  {
    micro->updateIntVars ();
  }

  void micropp_write_convergence_file_(void)
  {
    micro->writeConvergenceFile ();
  }
} 
