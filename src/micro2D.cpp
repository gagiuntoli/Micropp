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
template<>
int micropp<2>::get_elem_type(int ex, int ey)
{
	assert(micro_type == 0 || micro_type == 1);

	if (micro_type == 0) {
		// esfera en matriz
		const double x1 = ex * dx + dx / 2;
		const double y1 = ey * dy + dy / 2;
		const double x2 = lx / 2;
		const double y2 = ly / 2;
		const double rad = micro_params[3];
		return ((x2 - x1) * (x2 - x1) +
		        (y2 - y1) * (y2 - y1) < width * width);

	} else if (micro_type == 1) {
		const double y = ey * dy + dy / 2;
		return (y < width);
	}

	cerr << "Invalid micro_type = " << micro_type << endl;
	return -1;
}


template<>
micropp<2>::micropp(const int _ngp, const int size[3], const int _micro_type,
                    const double *_micro_params, const int *_mat_types,
                    const double *_params):
	ngp(_ngp),
	nx(size[0]), ny(size[1]), nz(1),
	nn(nx * ny),
	nex(nx - 1), ney(ny - 1), nez(1),
	nelem(nex * ney * nez),
	lx(_micro_params[0]), ly(_micro_params[1]),	lz(0.0),
	dx(lx / nex), dy(ly / ney), dz(0),
	dxi(1 / dx), dyi(1 / dy), dzi(0),

	width(_micro_params[3]), inv_tol(_micro_params[4]),
	micro_type(_micro_type), num_int_vars(nelem * 8 * NUM_VAR_GP)
{
	assert(dim == 2);

	initialize(_micro_params, _mat_types, _params);

	for (int ex = 0; ex < nx - 1; ex++)
		for (int ey = 0; ey < ny - 1; ey++)
			elem_type[glo_elem(ex, ey, 0)] = get_elem_type(ex, ey);

	const int ns[3] = { nx, ny, nz };
	const int nfield = dim;
	ell_init(&A, nfield, dim, ns, CG_MIN_ERR, CG_MAX_ITS);

	calc_ctan_lin();
}


template <>
void micropp<2>::set_displ_bc(const double *eps)
{
	const double eps_t[2][2] = {
		{ eps[0], 0.5 * eps[2] },
		{ 0.5 * eps[2], eps[1] } };

	// y = 0
	for (int i = 0; i < nx; i++) {
		const int n = nod_index(i, 0, 0);
		const double coor[2] = { i * dx, 0 };
		mvp_2(eps_t, coor, &u[n * dim]);
	}

	// y = ly
	for (int i = 0; i < nx; i++) {
		const int n = nod_index(i, ny - 1, 0);
		const double coor[2] = { i * dx, ly };
		mvp_2(eps_t, coor, &u[n * dim]);
	}

	// x = 0
	for (int j = 1; j < ny - 1; j++) {
		const int n = nod_index(0, j, 0);
		const double coor[2] = { 0, j * dy };
		mvp_2(eps_t, coor, &u[n * dim]);
	}

	// x = lx
	for (int j = 1; j < ny - 1; j++) {
		const int n = nod_index(nx - 1, j, 0);
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
template <>
void micropp<2>::calc_bmat(int gp, double bmat[6][3 * 8]) const
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
template <>
void micropp<2>::get_elem_rhs(double *be, int ex, int ey) const
{
	double bmat[6][3 * 8];
	const double wg = 0.25 * dx * dy;

	memset(be, 0, 2 * 4 * sizeof(double));

	for (int gp = 0; gp < 4; gp++) {

		calc_bmat(gp, bmat);

		double strain_gp[3], stress_gp[3];
		get_strain(gp, strain_gp, ex, ey);
		get_stress(gp, strain_gp, stress_gp, ex, ey);

		for (int i = 0; i < npe * dim; i++)
			for (int j = 0; j < nvoi; j++)
				be[i] += bmat[j][i] * stress_gp[j] * wg;

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
void micropp<2>::get_elem_mat(double Ae[2 * 4 * 2 * 4], int ex, int ey) const
{
	INST_START;

	const int e = glo_elem(ex, ey, 0);
	const material_t material = get_material(e);
	const double wg = (1 / 4.0) * dx * dy;

	const double E = material.E;
	const double nu = material.nu;
	const bool plasticity = material.plasticity;
	const int npedim = npe * dim;

	double ctan[3][3] = { { (1 - nu),     	nu,                0 },
	                      {       nu, (1 - nu),                0 },
	                      {        0,        0, (1 - 2 * nu) / 2 } };

	for (int i = 0; i < nvoi; i++)
		for (int j = 0; j < nvoi; j++)
			ctan[i][j] *= E / ((1 + nu) * (1 - 2 * nu));

	memset(Ae, 0, npedim * npedim * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		double bmat[6][3 * 8], cxb[3][2 * 4] = { 0.0 };
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
					Ae[inpedim + j] += bmatmi * cxb[m][j] * wg;
			}
		}
	}
}


template <>
void micropp<2>::assembly_mat()
{
	INST_START;

	ell_set_zero_mat(&A);

	double Ae[2 * 4 * 2 * 4];
	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			get_elem_mat(Ae, ex, ey);
			ell_add_2D(&A, ex, ey, Ae);
		}
	}
	ell_set_bc_2D(&A);
}


template <>
void micropp<2>::calc_ave_stress(double stress_ave[6]) const
{
	const double wg = 0.25 * dx * dy;
	memset(stress_ave, 0, nvoi * sizeof(double));

	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {

			double stress_aux[3] = { 0.0 };

			for (int gp = 0; gp < npe; ++gp) {

				double strain_gp[3], stress_gp[3];
				get_strain(gp, strain_gp, ex, ey);
				get_stress(gp, strain_gp, stress_gp, ex, ey);
				for (int v = 0; v < nvoi; ++v)
					stress_aux[v] += stress_gp[v] * wg;

			}
			for (int v = 0; v < nvoi; ++v)
				stress_ave[v] += stress_aux[v];
		}
	}

	for (int v = 0; v < nvoi; ++v)
		stress_ave[v] /= (lx * ly);
}


template <>
void micropp<2>::calc_ave_strain(double strain_ave[6]) const
{
	const double wg = 0.25 * dx * dy;
	memset(strain_ave, 0, nvoi * sizeof(double));

	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {

			double strain_aux[3] = { 0.0 };

			for (int gp = 0; gp < npe; ++gp) {

				double strain_gp[6];
				get_strain(gp, strain_gp, ex, ey);
				for (int v = 0; v < nvoi; ++v)
					strain_aux[v] += strain_gp[v] * wg;
			}

			for (int v = 0; v < nvoi; v++)
				strain_ave[v] += strain_aux[v];
		}
	}

	for (int v = 0; v < nvoi; v++)
		strain_ave[v] /= (lx * ly);
}


template <>
void micropp<2>::calc_fields()
{
	const double ivol = 1.0 / (dx * dy);
	const double wg = 0.25 * dx * dy;

	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {

			double strain_aux[3] = { 0.0 };
			double stress_aux[3] = { 0.0 };

			for (int gp = 0; gp < npe; ++gp) {

				double stress_gp[3], strain_gp[3];
				get_strain(gp, strain_gp, ex, ey);
				get_stress(gp, strain_gp, stress_gp, ex, ey);
				for (int v = 0; v < nvoi; v++) {
					strain_aux[v] += strain_gp[v] * wg;
					stress_aux[v] += stress_gp[v] * wg;
				}

			}

			const int e = glo_elem(ex, ey, 0);
			for (int v = 0; v < nvoi; v++) {
				elem_strain[e * nvoi + v] = strain_aux[v] * ivol;
				elem_stress[e * nvoi + v] = stress_aux[v] * ivol;
			}
		}
	}
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
