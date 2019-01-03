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

#include <cmath>
#include <cassert>

#include "instrument.hpp"
#include "micro.hpp"


template <int tdim>
void micropp<tdim>::set_macro_strain(const int gp_id,
				     const double *macro_strain)
{
	assert(gp_id >= 0);
	assert(gp_id < ngp);
	memcpy(gp_list[gp_id].macro_strain, macro_strain, nvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_macro_stress(const int gp_id,
				     double *macro_stress) const
{
	assert(gp_id >= 0);
	assert(gp_id < ngp);
	memcpy(macro_stress, gp_list[gp_id].macro_stress, nvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_macro_ctan(const int gp_id, double *macro_ctan) const
{
	assert(gp_id >= 0);
	assert(gp_id < ngp);
	memcpy(macro_ctan, gp_list[gp_id].macro_ctan, nvoi * nvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::homogenize()
{
	INST_START;

#pragma omp parallel for
	for (int igp = 0; igp < ngp; ++igp) {

		newton_t newton;
		gp_t<tdim> * const gp_ptr = &gp_list[igp];

		int thread_id = omp_get_thread_num();
		ell_matrix *A_ptr = &matrices_A[thread_id];
		ell_matrix *A0_ptr = &matrices_A0[thread_id];

		gp_ptr->sigma_cost = 0;

		double *vars_new;
		if (gp_ptr->allocated) {
			vars_new = gp_ptr->int_vars_k;
		} else {
			vars_new = vars_new_aux;
		}

		// SIGMA 1 Newton-Raphson
		memcpy(gp_ptr->u_k, gp_ptr->u_n, nndim * sizeof(double));

		newton_raphson_v(gp_ptr->allocated,
				 NR_MAX_ITS,
				 MAT_MODE_A,
				 A_ptr,
				 A0_ptr,
				 gp_ptr->macro_strain,
				 gp_ptr->int_vars_n,
				 gp_ptr->u_k,
				 &newton,
				 false);

		memcpy(&(gp_ptr->newton), &newton, sizeof(newton_t));

		for (int i = 0; i < newton.its; ++i)
			gp_ptr->sigma_cost += newton.solver_its[i];

		calc_ave_stress(gp_ptr->u_k, gp_ptr->int_vars_n, gp_ptr->macro_stress);
		filter(gp_ptr->macro_stress, nvoi, FILTER_REL_TOL);

		/* Updates <vars_new> and <f_trial_max> */
		bool nl_flag = calc_vars_new(gp_ptr->u_k, gp_ptr->int_vars_n, vars_new, &f_trial_max);

		if (nl_flag == true) {
			if (gp_ptr->allocated == false) {
				gp_ptr->allocate(num_int_vars);
				memcpy(gp_ptr->int_vars_k, vars_new, num_int_vars * sizeof(double));
			}
		}

		if (gp_ptr->allocated) {

			// CTAN 3/6 Newton-Raphsons in 2D/3D
			double eps_1[6], sig_0[6], sig_1[6];

			memcpy(u_aux, gp_ptr->u_k, nndim * sizeof(double));
			memcpy(sig_0, gp_ptr->macro_stress, nvoi * sizeof(double));

			for (int i = 0; i < nvoi; ++i) {

				memcpy(eps_1, gp_ptr->macro_strain, nvoi * sizeof(double));
				eps_1[i] += D_EPS_CTAN_AVE;

				newton_raphson_v(true,
						 NR_MAX_ITS,
						 MAT_MODE_A,
						 A_ptr,
						 A0_ptr,
						 eps_1,
						 gp_ptr->int_vars_n,
						 u_aux,
						 nullptr,
						 false);

				for (int i = 0; i < newton.its; ++i)
					gp_ptr->sigma_cost += newton.solver_its[i];

				calc_ave_stress(u_aux, gp_ptr->int_vars_n, sig_1);

				for (int v = 0; v < nvoi; ++v)
					gp_ptr->macro_ctan[v * nvoi + i] =
						(sig_1[v] - sig_0[v]) / D_EPS_CTAN_AVE;

			}
			filter(gp_ptr->macro_ctan, nvoi * nvoi, FILTER_REL_TOL);
		}
	}
}


template <int tdim>
void micropp<tdim>::update_vars()
{
	INST_START;

	for (int igp = 0; igp < ngp; ++igp)
		gp_list[igp].update_vars();
}


template class micropp<2>;
template class micropp<3>;
