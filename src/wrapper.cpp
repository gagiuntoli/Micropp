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

  void micropp_output_(int *time_step, int *Gauss_ID, double *MacroStrain)
  {
    micro->output (*time_step, *Gauss_ID, MacroStrain);
  }

  void micropp_get_non_linear_flag_(int *Gauss_ID, int *non_linear)
  {
    micro->getNonLinearFlag (*Gauss_ID, non_linear);
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

  void micropp_write_convergence_file_(void)
  {
    micro->writeConvergenceFile ();
  }
} 
