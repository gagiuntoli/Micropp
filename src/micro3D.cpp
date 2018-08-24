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


template <>
void micropp<3>::set_displ_bc(const double eps[nvoi], double *u)
{
	const double eps_t[dim][dim] = {
		{       eps[0], 0.5 * eps[3], 0.5 * eps[4] },
		{ 0.5 * eps[3],       eps[1], 0.5 * eps[5] },
		{ 0.5 * eps[4], 0.5 * eps[5],       eps[2] } };


	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, 0); // z = 0
			const double coor[3] = { i * dx, j * dy, 0 };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, nz - 1); // z = lz
			const double coor[3] = { i * dx, j * dy, lz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, 0, k); // y = 0
			const double coor[3] = { i * dx, 0, k * dz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, ny - 1, k); // y = ly
			const double coor[3] = { i * dx, ly , k * dz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int j = 1; j < ny - 1; ++j) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(0, j, k); // x = 0
			const double coor[3] = { 0, j * dy , k * dz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int j = 1; j < ny - 1; ++j) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(nx - 1, j, k); // x = lx
			const double coor[3] = { ly, j * dy , k * dz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}
}


template <>
void micropp<3>::calc_bmat(int gp, double bmat[nvoi][npe * dim]) const
{
	INST_START;

	const double dsh[8][3] = {
		{ -(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
		  -(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
		  -(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{ +(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
		  -(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
		  -(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{ +(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
		  +(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
		  -(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{ -(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
		  +(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
		  -(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{ -(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
		  -(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
		  +(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{ +(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
		  -(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
		  +(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{ +(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
		  +(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
		  +(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{ -(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
		  +(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
		  +(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz } };

	for (int i = 0; i < 8; ++i) {
		bmat[0][i * dim    ] = dsh[i][0];
		bmat[0][i * dim + 1] = 0;
		bmat[0][i * dim + 2] = 0;
		bmat[1][i * dim    ] = 0;
		bmat[1][i * dim + 1] = dsh[i][1];
		bmat[1][i * dim + 2] = 0;
		bmat[2][i * dim    ] = 0;
		bmat[2][i * dim + 1] = 0;
		bmat[2][i * dim + 2] = dsh[i][2];
		bmat[3][i * dim    ] = dsh[i][1];
		bmat[3][i * dim + 1] = dsh[i][0];
		bmat[3][i * dim + 2] = 0;
		bmat[4][i * dim    ] = dsh[i][2];
		bmat[4][i * dim + 1] = 0;
		bmat[4][i * dim + 2] = dsh[i][0];
		bmat[5][i * dim    ] = 0;
		bmat[5][i * dim + 1] = dsh[i][2];
		bmat[5][i * dim + 2] = dsh[i][1];
	}
}


template<>
double micropp<3>::assembly_rhs(const double *u)
{
	INST_START;

	memset(b, 0, nndim * sizeof(double));

	double be[dim * npe];
	int index[dim * npe];

	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				int n[npe];
				get_elem_nodes(n, ex, ey, ez);

				for (int j = 0; j < npe; ++j)
					for (int d = 0; d < dim; ++d)
						index[j * dim + d] = n[j] * dim + d;

				get_elem_rhs(u, be, ex, ey, ez);

				for (int i = 0; i < npe * dim; ++i)
					b[index[i]] += be[i];
			}
		}
	}

	// boundary conditions
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, 0); // z = 0
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, nz - 1); // z = lx
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, 0, k); // y = 0
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, ny - 1, k); // y = ly
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int j = 1; j < ny - 1; ++j) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(0, j, k); // x = 0
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int j = 1; j < ny - 1; j++) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(nx - 1, j, k); // x = lx
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	// Common
	for (int i = 0; i < nndim; ++i)
		b[i] = -b[i];

	double norm = 0.0;
	for (int i = 0; i < nndim; ++i)
		norm += b[i] * b[i];
	norm = sqrt(norm);

	return norm;
}


template <>
template <>
void micropp<3>::get_elem_mat(const double *u,
							  double Ae[npe * dim * npe * dim],
							  int ex, int ey, int ez) const
{
	INST_START;
	const int e = glo_elem(ex, ey, ez);
	const material_t material = get_material(e);

	double ctan[nvoi][nvoi];
	constexpr int npedim = npe * dim;
	constexpr int npedim2 = npedim * npedim;

	double TAe[npedim2] = { 0.0 };

	// memset(Ae, 0, npedim2 * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		double eps[6];
		get_strain(u, gp, eps, ex, ey, ez);
		const double *eps_p_old = &vars_old[intvar_ix(e, gp, 0)];
		double alpha_old = vars_old[intvar_ix(e, gp, 6)];

		if (material.plasticity)
			plastic_get_ctan(&material, eps, eps_p_old, alpha_old, ctan);
		else
			isolin_get_ctan(&material, ctan);

		double bmat[nvoi][npedim], cxb[nvoi][npedim];
		calc_bmat(gp, bmat);

		for (int i = 0; i < nvoi; ++i) {
			for (int j = 0; j < npedim; ++j) {
				double tmp = 0.0;
				for (int k = 0; k < nvoi; ++k)
					tmp += ctan[i][k] * bmat[k][j];
				cxb[i][j] = tmp;
			}
		}

		for (int m = 0; m < nvoi; ++m) {
			for (int i = 0; i < npedim; ++i) {
				const int inpedim = i * npedim;
				const double bmatmi = bmat[m][i];
				for (int j = 0; j < npedim; ++j)
					TAe[inpedim + j] += bmatmi * cxb[m][j] * wg;
			}
		}
		memcpy(Ae, TAe, npedim2 * sizeof(double));
	}
}


template <>
void micropp<3>::assembly_mat(const double *u)
{
	INST_START;

	ell_set_zero_mat(&A);

	double Ae[npe * dim * npe * dim];
	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {
				get_elem_mat(u, Ae, ex, ey, ez);
				ell_add_3D(&A, ex, ey, ez, Ae);
			}
		}
	}
	ell_set_bc_3D(&A);
}


// Explicit instantiation
template class micropp<3>;
