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
int micropp<tdim>::newton_raphson(double strain[nvoi],
				  double *int_vars_old,
				  double *u,
				  double newton_err[NR_MAX_ITS],
				  int solver_its[NR_MAX_ITS],
				  double solver_err[NR_MAX_ITS])
{
	INST_START;

	set_displ_bc(strain, u);

	int lits = 0;
	double lerr = 0.0, cg_err;

	if (solver_its != NULL) memset(solver_its, 0, NR_MAX_ITS * sizeof(int));
	if (solver_err != NULL) memset(solver_err, 0, NR_MAX_ITS * sizeof(double));
	if (newton_err != NULL) memset(newton_err, 0, NR_MAX_ITS * sizeof(double));

	do {
		lerr = assembly_rhs(u, int_vars_old);

		if (newton_err != NULL) newton_err[lits] = lerr;

		if (lerr < NR_MAX_TOL)
			break;

		assembly_mat(&A, u, int_vars_old);

		int cg_its = ell_solve_cgpd(&A, b, du, &cg_err);

		if (solver_its != NULL) solver_its[lits] = cg_its;
		if (solver_err != NULL) solver_err[lits] = cg_err;

		for (int i = 0; i < nn * dim; ++i)
			u[i] += du[i];

		lits++;

	} while (lits < NR_MAX_ITS);

	return lits;
}


template class micropp<2>;
template class micropp<3>;
