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


template<int tdim>
void micropp<tdim>::homogenize_task_acc(int igp)
{
	const int ns[3] = { nx, ny, nz };

	ell_matrix A;  // Jacobian
	ell_init(&A, dim, dim, ns, CG_MIN_ERR, CG_REL_ERR, CG_MAX_ITS);
	double *b = (double *) calloc(nndim, sizeof(double));
	double *du = (double *) calloc(nndim, sizeof(double));
	double *u = (double *) calloc(nndim, sizeof(double));
	double *vars_new_aux = (double *) calloc(nvars, sizeof(double));

	gp_t<tdim> * const gp_ptr = &gp_list[igp];

	double *vars_new = (gp_ptr->allocated) ? gp_ptr->vars_k : vars_new_aux;

	// SIGMA 1 Newton-Raphson
	memcpy(u, gp_ptr->u_n, nndim * sizeof(double));

	newton_t newton = newton_raphson_acc(&A, b, u, du, gp_ptr->strain, gp_ptr->vars_n);

	memcpy(gp_ptr->u_k, u, nndim * sizeof(double));
	gp_ptr->cost += newton.solver_its;
	gp_ptr->converged = newton.converged;

	/*
	 * In case it has not converged do the sub-iterations
	 */
	if (gp_ptr->converged == false && subiterations == true) {

		double eps_sub[nvoi], deps_sub[nvoi];
		memcpy(u, gp_ptr->u_n, nndim * sizeof(double));
		memcpy(eps_sub, gp_ptr->strain_old, nvoi * sizeof(double));
		gp_ptr->subiterated = true;

		for (int i = 0; i < nvoi; ++i)
			deps_sub[i] = (gp_ptr->strain[i] - gp_ptr->strain_old[i]) / nsubiterations;

		for (int its = 0; its < nsubiterations; ++its) {

			for (int j = 0; j < nvoi; ++j)
				eps_sub[j] += deps_sub[j];

			newton = newton_raphson_acc(&A, b, u, du, eps_sub, gp_ptr->vars_n);
			gp_ptr->cost += newton.solver_its;
		}

		gp_ptr->converged = newton.converged;
		memcpy(gp_ptr->u_k, u, nndim * sizeof(double));
	}

	if (coupling == ONE_WAY) {

		memset (gp_ptr->stress, 0.0, nvoi * sizeof(double));
		for (int i = 0; i < nvoi; ++i)
			for (int j = 0; j < nvoi; ++j)
				gp_ptr->stress[i] += ctan_lin[i * nvoi + j] * gp_ptr->strain[j];

	} else if (coupling == FULL || coupling == NO_COUPLING) {

		calc_ave_stress(gp_ptr->u_k, gp_ptr->stress, gp_ptr->vars_n);
		filter(gp_ptr->stress, nvoi, FILTER_REL_TOL);

	}

	// Updates <vars_new>
	bool non_linear = calc_vars_new_acc(gp_ptr->u_k, gp_ptr->vars_n, vars_new);

	if (non_linear == true) {
		if (gp_ptr->allocated == false) {
			gp_ptr->allocate(nvars);
			memcpy(gp_ptr->vars_k, vars_new, nvars * sizeof(double));
		}
	}

	if (gp_ptr->allocated && coupling == FULL) {

		// CTAN 3/6 Newton-Raphsons in 2D/3D
		double eps_1[6], sig_0[6], sig_1[6];

		memcpy(u, gp_ptr->u_k, nndim * sizeof(double));
		memcpy(sig_0, gp_ptr->stress, nvoi * sizeof(double));

		for (int i = 0; i < nvoi; ++i) {

			memcpy(eps_1, gp_ptr->strain, nvoi * sizeof(double));
			eps_1[i] += D_EPS_CTAN_AVE;

			newton_raphson_acc(&A, b, u, du, eps_1, gp_ptr->vars_n);

			gp_ptr->cost += newton.solver_its;

			calc_ave_stress(u, sig_1, gp_ptr->vars_n);

			for (int v = 0; v < nvoi; ++v)
				gp_ptr->ctan[v * nvoi + i] =
					(sig_1[v] - sig_0[v]) / D_EPS_CTAN_AVE;

		}
		filter(gp_ptr->ctan, nvoi * nvoi, FILTER_REL_TOL);
	}

	ell_free(&A);
	free(b);
	free(u);
	free(du);
	free(vars_new_aux);
}


template class micropp<3>;
