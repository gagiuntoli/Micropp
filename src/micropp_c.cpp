/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
 *                         Judicaël Grasset <judicael.grasset@stfc.ac.uk>
 *                         Alejandro Figueroa <afiguer7@maisonlive.gmu.edu>
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
			  const int type, const double *geo_params,
	                  const struct material_base *materials,
			  const int *coupling, const int nsubiterations,
			  const int mpi_rank)
	{
		micropp_params_t params;
		params.ngp = ngp;
		memcpy(params.size, size, 3 * sizeof(int));
		params.type = type;
		memcpy(params.geo_params, geo_params, 4 * sizeof(double));
		for (int i = 0; i < MAX_MATERIALS; ++i) {
			params.materials[i] = materials[i];
		}
		params.coupling = new int[ngp];
		memcpy(params.coupling, coupling, ngp * sizeof(int));
		params.subiterations = true;
		params.nsubiterations = nsubiterations;
		params.mpi_rank = mpi_rank;
		params.use_A0 = false;
		params.its_with_A0 = 1;

		self->ptr = new micropp<3>(params);
	}

	void micropp3_free(micropp3 *self)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		delete ptr;
	}

	void micropp3_set_strain(micropp3 *self, const int gp_id,
				 const double *strain)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->set_strain(gp_id, strain);
	}

	void micropp3_get_stress(const micropp3 *self, const int gp_id,
				 double *stress)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->get_stress(gp_id, stress);
	}

	void micropp3_get_ctan(const micropp3 *self, int gp,
			       double ctan[36])
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->get_ctan(gp, ctan);
	}

	void micropp3_homogenize(micropp3 *self, const int homog_type)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->homogenize(homog_type);
	}

	int micropp3_get_cost(const micropp3 *self, int gp_id)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		return ptr->get_cost(gp_id);
	}

	bool micropp3_has_converged(const micropp3 *self, const int gp_id)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		return ptr->has_converged(gp_id);
	}

	bool micropp3_has_subiterated(const micropp3 *self, const int gp_id)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		return ptr->has_subiterated(gp_id);
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

	void micropp3_output2(micropp3 *self, const int gp_id,
			      const int elem_global, const int time_step)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->output2(gp_id, elem_global, time_step);
	}

	void micropp3_print_info(micropp3 *self)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->print_info();
	}

	bool micropp3_is_non_linear(const micropp3 *self, const int gp_id)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		return ptr->is_non_linear(gp_id);
	}

	int micropp3_get_non_linear_gps(const micropp3 *self)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		return ptr->get_non_linear_gps();
	}

	void micropp3_write_restart(const micropp3 *self, const int restart_id)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->write_restart(restart_id);
	}

	void micropp3_read_restart(const micropp3 *self, const int restart_id)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->read_restart(restart_id);
	}

	void micropp3_write_profiling(const micropp3 *self, const int profile_id)
	{
		micropp<3> *ptr = (micropp<3> *) self->ptr;
		ptr->write_profiling(profile_id);
	}
}
