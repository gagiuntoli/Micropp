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

#include "micropp_c.h"
#include "micro.hpp"
#include "material_base.h"

extern "C" {

	// IMPORTANT!! This struct should match with the one in FORTRAN

	void micropp3_new(struct micropp3 *self, int ngp, const int size[3],
	                  const int micro_type, const double *micro_params,
	                  const material_base *materials)
	{
		material_t *tmp = (material_t *) materials;
		self->ptr = new micropp<3>(ngp, size, micro_type,
		                           micro_params, tmp);
	}

	void micropp3_free(micropp3 *self)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		delete ptr;
	}

	void micropp3_set_macro_strain(micropp3 *self, const int gp_id,
				const double *macro_strain)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->set_macro_strain(gp_id, macro_strain);
	}

	void micropp3_homogenize(micropp3 *self)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->homogenize();
	}

	int micropp3_get_cost(const micropp3 *self, int gp_id)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		return ptr->get_cost(gp_id);
	}

	void micropp3_get_macro_stress(const micropp3 *self, const int gp_id,
	                                 double *macro_stress)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->get_macro_stress(gp_id, macro_stress);
	}

	void micropp3_update_vars(micropp3 *self)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->update_vars();
	}

	void micropp3_output(micropp3 *self, const int gp_id,
	                        const char *filename)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->output(gp_id, filename);
	}

	void micropp3_print_info(micropp3 *self)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->print_info();
	}

	bool micropp3_get_nl_flag(const micropp3 *self, const int gp_id)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		return ptr->is_non_linear(gp_id);
	}

	int micropp3_get_non_linear_gps(const micropp3 *self)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		return ptr->get_non_linear_gps();
	}

	double micropp3_get_f_trial_max(const micropp3 *self)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		return ptr->get_f_trial_max();
	}

	void micropp3_get_macro_ctan(const micropp3 *self, int gp, double ctan[36])
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->get_macro_ctan(gp, ctan);
	}



}
