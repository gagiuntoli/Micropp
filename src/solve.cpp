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
int micropp<tdim>::newton_raphson_linear(ell_matrix *A0,
					 double *b,
					 double *u,
					 double *du,
					 const double strain[nvoi],
					 bool print)






{

	int thread_id = omp_get_thread_num();
	return newton_raphson_v(nullptr,
				A0,
				b,
				u,
				du,
				false,
				2,
				MAT_MODE_A0,
				strain,
				nullptr,
				nullptr,
				print);
}


template <int tdim>
int micropp<tdim>::newton_raphson_v(ell_matrix *A,
				    ell_matrix *A0,
				    double *b,
				    double *u,
				    double *du,
				    const bool non_linear,
				    const int newton_max_its,
				    const int mat_mode,
				    const double strain[nvoi],
				    const double *vars_old,
				    newton_t *newton,
				    bool print)
{
	INST_START;

	set_displ_bc(strain, u);

	int its = 0;
	double norm, norm_0, cg_err;
	ell_matrix *A_ptr;

	if (newton != nullptr) {
		memset(newton->solver_its, 0, NR_MAX_ITS * sizeof(int));
		memset(newton->solver_norms, 0, NR_MAX_ITS * sizeof(double));
		memset(newton->norms, 0, NR_MAX_ITS * sizeof(double));
	}

	int thread_id = omp_get_thread_num();

	while (its < newton_max_its) {

		norm = assembly_rhs(u, vars_old, b);
		if (print)
			printf("Thread : %d |b| = %e\n", thread_id, norm);

		if (its == 0)
			norm_0 = norm;

		if (newton != nullptr)
			newton->norms[its] = norm;

		if (norm < NR_MAX_TOL || norm < norm_0 * NR_REL_TOL)
			break;

		int cg_its;
		switch (mat_mode) {

			case MAT_MODE_A0:

				/* Always use A0 */
				A_ptr = A0;
				break;

			case MAT_MODE_A:

				/* Assemblies a new matrix if we are in a
				 * non-linear situation or if we are in a
				 * second newton-raphson iteration.
				 */

				if (non_linear || its > 0) {

					assembly_mat(A, u, vars_old);
					A_ptr = A;

				} else {

					A_ptr = A0;

				}

				break;

			default:
				return 1;
				break;
		}

		if (print)
			printf("Thread : %d SOLVER_START\n", thread_id);

		cg_its = ell_solve_cgpd(A_ptr, b, du, &cg_err);

		if (print)
			printf("Thread : %d SOLVER_END ITS : %d\n", thread_id, cg_its);

		if (newton != nullptr) {
			newton->solver_its[its] = cg_its;
			newton->solver_norms[its] = cg_err;
		}

		for (int i = 0; i < nn * dim; ++i)
			u[i] += du[i];

		its++;

	}

	if (newton != nullptr)
		newton->its = its;

	return 0;
}


template class micropp<2>;
template class micropp<3>;
