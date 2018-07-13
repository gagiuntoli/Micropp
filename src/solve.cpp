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

#include <iostream>
#include <iomanip>		// print with format

#include <cassert>

#include "instrument.hpp"
#include "micro.hpp"

using namespace std;

void micropp_t::solve()
{
	INST_START;

	assert(dim == 2 | dim ==3);
	ell_solve_cgpd_struct(&solver, &A, dim, dim, nn, b, du);
}

void micropp_t::newton_raphson(bool *nl_flag, int *its, double *err)
{
	INST_START;

	int lits = 0;
	double lerr = 0.0;

	do {
		lerr = assembly_rhs(nl_flag);

		if (lerr < NR_MAX_TOL)
			break;

		assembly_mat();

		memset(du, 0.0, nn * dim * sizeof(double));
		solve();

		for (int i = 0; i < nn * dim; ++i)
			u[i] += du[i];

		lits++;

	} while ((lits < NR_MAX_ITS) && (lerr > NR_MAX_TOL));

	*its = lits;
	*err = lerr;
}
