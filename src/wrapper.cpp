/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
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
#include "micro.hpp"

using namespace std;

static micropp_t* micro = NULL;

extern "C" {
	void micropp_construct_(int *dim, int size[3], int *micro_type,
	                        double *micro_params, int *mat_types, double *params)
	{
		micro = new micropp_t(*dim, size, *micro_type, micro_params, mat_types, params); 
	}

	void micropp_output_(int *tstep, int *gp_id)
	{
		micro->output (*tstep, *gp_id);
	}

	void micropp_get_non_linear_flag_(int *gp_id, int *non_linear)
	{
		micro->get_nl_flag (*gp_id, non_linear);
	}

	void micropp_set_macro_strain_(int *gp_id, double *macro_strain)
	{
		micro->set_macro_strain(*gp_id, macro_strain);
	}

	void micropp_homogenize_(void)
	{
		micro->homogenize();
	}

	void micropp_get_macro_stress_(int *gp_id, double *macro_stress)
	{
		micro->get_macro_stress(*gp_id, macro_stress);
	}

	void micropp_get_macro_ctan_(int *gp_id, double *macro_ctan)
	{
		micro->get_macro_ctan(*gp_id, macro_ctan);
	}

	void micropp_update_internal_variables_(void)
	{
		micro->update_vars ();
	}

	void micropp_write_convergence_file_(void)
	{
		micro->write_info_files ();
	}
}
