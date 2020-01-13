/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018
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
#include "common.hpp"


CUDA_HOSTDEV
#pragma acc routine seq
void get_elem_nodes(int n[8], const int nx, const int ny,
		    const int ex, const int ey, const int ez)
{
	const int nxny = ny * nx;
	const int n0 = ez * nxny + ey * nx + ex;
	n[0] = n0;
	n[1] = n0 + 1;
	n[2] = n0 + nx + 1;
	n[3] = n0 + nx;
	n[4] = n[0] + nxny;
	n[5] = n[1] + nxny;
	n[6] = n[2] + nxny;
	n[7] = n[3] + nxny;
}


CUDA_HOSTDEV
#pragma acc routine seq
template <int tdim>
void micropp<tdim>::get_elem_displ(const double *u,
				   double elem_disp[npe * dim],
				   int ex, int ey, int ez) const
{
	int n[npe];
	get_elem_nodes(n, nx, ny, ex, ey, ez);

	for (int i = 0 ; i < npe; ++i) {
		for (int d = 0; d < dim; ++d) {
			elem_disp[i * dim + d] = u[n[i] * dim + d];
		}
	}
}


CUDA_HOSTDEV
#pragma acc routine seq
template <int tdim>
void micropp<tdim>::get_strain(const double *u, int gp, double *strain_gp,
			       int ex, int ey, int ez) const
{
	double elem_disp[npe * dim];
	get_elem_displ(u, elem_disp, ex, ey, ez);

	for (int i = 0; i < nvoi; ++i) {
		strain_gp[i] = 0;
	}

	for (int v = 0; v < nvoi; ++v) {
		for (int i = 0; i < npe * dim; ++i){
			strain_gp[v] += bmat_cache[gp][v][i] * elem_disp[i];
		}
	}
}


template class micropp<3>;
