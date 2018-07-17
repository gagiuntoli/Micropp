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
void micropp<tdim>::newton_raphson(bool *nl_flag, int *its, double *err)
{
	INST_START;

	int lits = 0;
	double lerr = 0.0;

	do {
		lerr = assembly_rhs(nl_flag);  // Acts on b

		if (lerr < NR_MAX_TOL)
			break;

		assembly_mat();   // Acts on A

		// Move this inside ell_solve_cgpd_struct
		memset(du, 0.0, nn * dim * sizeof(double));

		// in(b) inout
		ell_solve_cgpd_struct(&solver, &A, dim, dim, nn, b, du);

		for (int i = 0; i < nn * dim; ++i)
			u[i] += du[i];

		lits++;

	} while (lits < NR_MAX_ITS);

	*its = lits;
	*err = lerr;
}


// Explicit instantiation
template class micropp<2>;
template class micropp<3>;
