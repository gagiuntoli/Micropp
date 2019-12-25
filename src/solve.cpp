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


#include "micropp.hpp"


using namespace std;


template <int tdim>
newton_t micropp<tdim>::newton_raphson(ell_matrix *A, double *b, double *u,
				       double *du, const double strain[nvoi],
				       const double *vars_old)
{

	INST_START;

	newton_t newton;

	set_displ_bc(strain, u);

	int its = 0;

#ifdef _OPENACC
	double norm = assembly_rhs_acc(u, vars_old, b);
#else
	double norm = assembly_rhs(u, vars_old, b);
#endif
	const double norm_0 = norm;

	while (its < nr_max_its) {

		if (norm < nr_max_tol || norm < norm_0 * nr_rel_tol) {
			newton.converged = true;
			break;
		}

		/*
		 * Matrix selection according if it's linear or non-linear.
		 * All OpenMP threads can access to A0 with no cost because
		 * is a read-only matrix.
		 *
		 */
		ell_matrix *A_ptr;
		if (!use_A0 || its > (its_with_A0 - 1)) {
#if defined (_OPENACC)
			assembly_mat_acc(A, u, vars_old);
			cout << "solving with OPENACC" << endl;
#elif defined (_CUDA)
			assembly_mat_cuda(A, u, vars_old);
			cout << "solving with CUDA" << endl;
#else
			assembly_mat(A, u, vars_old);
			cout << "solving with CPU" << endl;
#endif
			A_ptr = A;
		} else {
#ifdef _OPENMP
			int tid = omp_get_thread_num();
#else
			int tid = 0;
#endif
			A_ptr = &A0[tid];
		}

		double cg_err;
#if defined (_OPENACC)
		int cg_its = ell_solve_cgpd_acc(A_ptr, b, du, &cg_err);
#elif defined(_CUDA)
		int cg_its = ell_solve_cgpd_cuda(A_ptr, b, du, &cg_err);
#else
		int cg_its = ell_solve_cgpd(A_ptr, b, du, &cg_err);
#endif

		newton.solver_its += cg_its;

		for (int i = 0; i < nn * dim; ++i)
			u[i] += du[i];

#ifdef _OPENACC
		norm = assembly_rhs_acc(u, vars_old, b);
#else
		norm = assembly_rhs(u, vars_old, b);
#endif

		its++;
	}

	newton.its = its;
	return newton;
}


template class micropp<3>;
