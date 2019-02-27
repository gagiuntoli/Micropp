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


#include "micro.hpp"


using namespace std;


template <int tdim>
int micropp<tdim>::newton_raphson(ell_matrix *A, double *b, double *u, double *du,
				  const bool non_linear, const double strain[nvoi],
				  const double *vars_old, newton_t *newton)
{
	INST_START;

	if (newton == nullptr)
		return -1;

	set_displ_bc(strain, u);

	int its = 0;

	memset(newton->solver_its, 0, NR_MAX_ITS * sizeof(int));
	memset(newton->solver_norms, 0, NR_MAX_ITS * sizeof(double));
	memset(newton->norms, 0, NR_MAX_ITS * sizeof(double));

	int thread_id = omp_get_thread_num();
	double norm = assembly_rhs(u, vars_old, b);
	const double norm_0 = norm;

	while (its < newton->max_its) {

		newton->norms[its] = norm;

		if (norm < newton->max_tol || norm < norm_0 * newton->rel_tol)
			break;

		assembly_mat(A, u, vars_old);

		double cg_err;
		int cg_its = ell_solve_cgpd(A, b, du, &cg_err);

		newton->solver_its[its] = cg_its;
		newton->solver_norms[its] = cg_err;

		for (int i = 0; i < nn * dim; ++i)
			u[i] += du[i];

		norm = assembly_rhs(u, vars_old, b);

		its++;
	}

	if (newton != nullptr)
		newton->its = its;

	return 0;
}


template class micropp<2>;
template class micropp<3>;
