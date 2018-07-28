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


template<>
micropp<2>::micropp(const int _ngp, const int size[3], const int _micro_type,
                    const double *_micro_params, const material_t *_materials):
	ngp(_ngp),
	nx(size[0]), ny(size[1]), nz(1),
	nn(nx * ny),
	nex(nx - 1), ney(ny - 1), nez(1),
	nelem(nex * ney * nez),
	lx(_micro_params[0]), ly(_micro_params[1]),	lz(0.0),
	dx(lx / nex), dy(ly / ney), dz(0.0),

	width(_micro_params[3]), inv_tol(_micro_params[4]),
	wg(dx * dy / npe), ivol(1.0 / (dx * dy)),
	micro_type(_micro_type), num_int_vars(nelem * 8 * NUM_VAR_GP)
{
	initialize(_micro_params, _materials);
}


template <>
void micropp<2>::set_displ_bc(const double *eps)
{
	const double eps_t[dim][dim] =
		{ {       eps[0], 0.5 * eps[2] },
		  { 0.5 * eps[2],       eps[1] } };

	for (int i = 0; i < nx; ++i) {
		const int n = nod_index(i, 0, 0); // y = 0
		const double coor[2] = { i * dx, 0 };
		mvp_2(eps_t, coor, &u[n * dim]);
	}

	for (int i = 0; i < nx; ++i) {
		const int n = nod_index(i, ny - 1, 0); // y = ly
		const double coor[2] = { i * dx, ly };
		mvp_2(eps_t, coor, &u[n * dim]);
	}

	for (int j = 1; j < ny - 1; ++j) {
		const int n = nod_index(0, j, 0); // x = 0
		const double coor[2] = { 0, j * dy };
		mvp_2(eps_t, coor, &u[n * dim]);
	}

	for (int j = 1; j < ny - 1; ++j) {
		const int n = nod_index(nx - 1, j, 0); // x = lx
		const double coor[2] = { lx, j * dy };
		mvp_2(eps_t, coor, &u[n * dim]);
	}
}


template <>
void micropp<2>::isolin_get_stress(
		const material_t *material, const double eps[6],
		double stress[6]) const
{
	const double E = material->E;
	const double nu = material->nu;

	double ctan[3][3] = { { (1 - nu),     	nu,                0 },
	                      {       nu, (1 - nu),                0 },
	                      {        0,        0, (1 - 2 * nu) / 2 } };

	for (int i = 0; i < nvoi; i++)
		for (int j = 0; j < nvoi; j++)
			ctan[i][j] *= E / ((1 + nu) * (1 - 2 * nu));

	memset(stress, 0, nvoi * sizeof(double));
	for (int i = 0; i < nvoi; i++)
		for (int j = 0; j < nvoi; j++)
			stress[i] += ctan[i][j] * eps[j];
}


template <>
void micropp<2>::calc_bmat(int gp, double bmat[nvoi][npe *dim])	const
{
	const double dsh[4][2] = {
		{ -(1 - xg[gp][1]) / 4 * 2 / dx,
		  -(1 - xg[gp][0]) / 4 * 2 / dy },
		{ +(1 - xg[gp][1]) / 4 * 2 / dx,
		  -(1 + xg[gp][0]) / 4 * 2 / dy },
		{ +(1 + xg[gp][1]) / 4 * 2 / dx,
		  +(1 + xg[gp][0]) / 4 * 2 / dy },
		{ -(1 + xg[gp][1]) / 4 * 2 / dx,
		  +(1 - xg[gp][0]) / 4 * 2 / dy } };

	for (int i = 0; i < 4; ++i) {
		bmat[0][i * dim    ] = dsh[i][0];
		bmat[0][i * dim + 1] = 0;
		bmat[1][i * dim    ] = 0;
		bmat[1][i * dim + 1] = dsh[i][1];
		bmat[2][i * dim    ] = dsh[i][1];
		bmat[2][i * dim + 1] = dsh[i][0];
	}
}


template <>
double micropp<2>::assembly_rhs()
{
	INST_START;

	memset(b, 0.0, nn * dim * sizeof(double));

	double be[2 * 4];
	int index[2 * 4];

	for (int ex = 0; ex < nx - 1; ++ex) {
		for (int ey = 0; ey < ny - 1; ++ey) {

			int n[8];
			get_elem_nodes(n, ex, ey);

			for (int j = 0; j < npe; ++j)
				for (int d = 0; d < dim; ++d)
					index[j * dim + d] = n[j] * dim + d;

			get_elem_rhs(be, ex, ey);

			for (int i = 0; i < npe * dim; ++i)
				b[index[i]] += be[i];

		}
	}

	// boundary conditions
	for (int i = 0; i < nx; i++) {
		const int n = nod_index2D(i, 0); // y = 0
		memset(&b[n * dim], 0, dim * sizeof(double));
	}
	for (int i = 0; i < nx; i++) {
		const int n = nod_index2D(i, ny - 1); // y = ly
		memset(&b[n * dim], 0, dim * sizeof(double));
	}
	for (int j = 1; j < ny - 1; j++) {
		const int n = nod_index2D(0, j); // x = 0
		memset(&b[n * dim], 0, dim * sizeof(double));
	}
	for (int j = 1; j < ny - 1; j++) {
		const int n = nod_index2D(nx - 1, j); // x = lx
		memset(&b[n * dim], 0, dim * sizeof(double));
	}

	// Common part
	for (int i = 0; i < nn * dim; ++i)
		b[i] = -b[i];

	double norm = 0.0;
	for (int i = 0; i < nn * dim; ++i)
		norm += b[i] * b[i];
	norm = sqrt(norm);

	return norm;
}


template <>
template <>
void micropp<2>::get_elem_mat(double Ae[npe * dim * npe * dim],
                              int ex, int ey) const
{
	INST_START;
	const int e = glo_elem(ex, ey, 0);
	const material_t material = get_material(e);

	const double E = material.E;
	const double nu = material.nu;
	const bool plasticity = material.plasticity;
	constexpr int npedim = npe * dim;
	constexpr int npedim2 = npedim * npedim;

	double ctan[nvoi][nvoi] = {{(1 - nu),     nu,             0 },
	                           {      nu, (1 - nu),           0 },
	                           {       0,   0, (1 - 2 * nu) / 2 }};

	for (int i = 0; i < nvoi; i++)
		for (int j = 0; j < nvoi; j++)
			ctan[i][j] *= E / ((1 + nu) * (1 - 2 * nu));

	double TAe[npedim2] = { 0.0 };

	for (int gp = 0; gp < npe; ++gp) {

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
void micropp<2>::assembly_mat()
{
	INST_START;

	ell_set_zero_mat(&A);

	double Ae[npe * dim * npe * dim];
	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			get_elem_mat(Ae, ex, ey);
			ell_add_2D(&A, ex, ey, Ae);
		}
	}
	ell_set_bc_2D(&A);
}


template<>
bool micropp<2>::calc_vars_new()
{
    return false;
}


template <>
void micropp<2>::plastic_get_stress(
		const material_t *material, const double eps[6],
		const double eps_p_old[6], double alpha_old,
		double stress[6]) const
{
}

// Explicit instantiation
template class micropp<2>;
