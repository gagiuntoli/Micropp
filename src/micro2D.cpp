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
	dx(lx / nex), dy(ly / nex), dz(0),

	width(_micro_params[3]), inv_tol(micro_params[4]),
	micro_type(_micro_type), num_int_vars(nelem * 8 * NUM_VAR_GP)
{
	assert(dim == 2);

	initialize(_micro_params, _mat_types, _params);

	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {
			int e = glo_elem3D(ex, ey, 0);
			elem_type[e] = get_elem_type(ex, ey);
		}
	}
	ell_init_2D(&A, dim, nx, ny);

	calc_ctan_lin();
}


template <>
void micropp<2>::set_displ(double *eps)
{
	// y = 0
	for (int i = 0; i < nx; i++) {
		const double xcoor = i * dx;
		const double ycoor = 0.0;
		const double dux = eps[0] * xcoor + 0.5 * eps[2] * ycoor;
		const double duy = 0.5 * eps[2] * xcoor + eps[1] * ycoor;
		u[i * dim] = dux;
		u[i * dim + 1] = duy;
	}
	// y = ly
	for (int i = 0; i < nx; i++) {
		const double xcoor = i * dx;
		const double ycoor = ly;
		const double dux = eps[0] * xcoor + 0.5 * eps[2] * ycoor;
		const double duy = 0.5 * eps[2] * xcoor + eps[1] * ycoor;
		u[(i + (ny - 1) * nx) * dim] = dux;
		u[(i + (ny - 1) * nx) * dim + 1] = duy;
	}
	// x = 0
	for (int i = 0; i < ny - 2; i++) {
		const double xcoor = 0.0;
		const double ycoor = (i + 1) * dy;
		const double dux = eps[0] * xcoor + 0.5 * eps[2] * ycoor;
		const double duy = 0.5 * eps[2] * xcoor + eps[1] * ycoor;
		u[(i + 1) * nx * dim] = dux;
		u[(i + 1) * nx * dim + 1] = duy;
	}
	// x = lx
	for (int i = 0; i < ny - 2; i++) {
		const double xcoor = lx;
		const double ycoor = (i + 1) * dy;
		const double dux = eps[0] * xcoor + 0.5 * eps[2] * ycoor;
		const double duy = 0.5 * eps[2] * xcoor + eps[1] * ycoor;
		u[((i + 2) * nx - 1) * dim] = dux;
		u[((i + 2) * nx - 1) * dim + 1] = duy;
	}
}


template <>
template <>
void micropp<2>::get_stress(int gp, double strain_gp[3], bool * non_linear,
                            double *stress_gp, int ex, int ey) const
{
	*non_linear = false;

	const int e = glo_elem3D(ex, ey, 0);
	const material_t material = get_material(e);

	const double E = material.E;
	const double nu = material.nu;

	double ctan[3][3] = { { (1 - nu),     	nu,                0 },
	                      {       nu, (1 - nu),                0 },
	                      {        0,        0, (1 - 2 * nu) / 2 } };

	for (int i = 0; i < nvoi; i++)
		for (int j = 0; j < nvoi; j++)
			ctan[i][j] *= E / ((1 + nu) * (1 - 2 * nu));

	for (int i = 0; i < 3; i++) {
		stress_gp[i] = 0.0;
		for (int j = 0; j < 3; j++) {
			stress_gp[i] += ctan[i][j] * strain_gp[j];
		}
	}
}


template <>
template <>
void micropp<2>::calc_bmat(int gp, double bmat[3][2 * 4]) const
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
void micropp<2>::getElemDisp(double *elem_disp, int ex, int ey) const
{
	const int n[] = { ey * nx + ex,
	                  ey * nx + ex + 1,
	                  (ey + 1) * nx + ex + 1,
	                  (ey + 1) * nx + ex };

	for (int i = 0 ; i < npe; ++i)
		for (int d = 0; d < dim; ++d)
			elem_disp[i * dim + d] = u[n[i] * dim + d];

}


template <>
template <>
void micropp<2>::get_strain(int gp, double *strain_gp,
                            int ex, int ey) const
{
	double elem_disp[2 * 4];
	getElemDisp(elem_disp, ex, ey);

	double bmat[3][8];
	calc_bmat(gp, bmat);

	for (int v = 0; v < nvoi; ++v) {
		strain_gp[v] = 0.0;
		for (int i = 0; i < npe * dim; i++) {
			strain_gp[v] += bmat[v][i] * elem_disp[i];
		}
	}
}


template <>
template <>
void micropp<2>::get_elem_rhs(bool *non_linear,
                              double *be, int ex, int ey) const
{
	double bmat[3][2 * 4], cxb[3][8], stress_gp[6];
	const double wg = 0.25 * dx * dy;

	for (int i = 0; i < 2 * 4; i++)
		be[i] = 0.0;

	for (int gp = 0; gp < 4; gp++) {

		calc_bmat(gp, bmat);

		double strain_gp[3];
		get_strain(gp, strain_gp, ex, ey);
		get_stress(gp, strain_gp, non_linear, stress_gp, ex, ey);

		for (int i = 0; i < npe * dim; i++) {
			for (int j = 0; j < nvoi; j++) {
				be[i] += bmat[j][i] * stress_gp[j] * wg;
			}
		}

	}			// gp loop
}


template <>
double micropp<2>::assembly_rhs(bool *non_linear)
{
	INST_START;

	bool non_linear_one_elem = false;
	*non_linear = false;

	for (int i = 0; i < nn * dim; ++i)
		b[i] = 0.0;

	double be[2 * 4];
	int index[2 * 4];

	for (int ex = 0; ex < nx - 1; ++ex) {
		for (int ey = 0; ey < ny - 1; ++ey) {

			const int n[npe] = { ey * nx + ex,
			                   ey * nx + ex + 1,
			                   (ey + 1) * nx + ex + 1,
			                   (ey + 1) * nx + ex };

			for (int j = 0; j < npe; ++j) {
				for (int d = 0; d < dim; ++d) {
					index[j * dim + d] = n[j] * dim + d;
				}
			}

			get_elem_rhs(&non_linear_one_elem, be, ex, ey);
			if (non_linear_one_elem)
				*non_linear = true;

			for (int i = 0; i < npe * dim; ++i)
				b[index[i]] += be[i];	// assembly

		}
	}

	// boundary conditions
	// y = 0
	for (int i = 0; i < nx; i++)
		for (int d = 0; d < dim; d++)
			b[i * dim + d] = 0.0;

	// y = ly
	for (int i = 0; i < nx; i++)
		for (int d = 0; d < dim; d++)
			b[(i + (ny - 1) * nx) * dim + d] = 0.0;

	// x = 0
	for (int i = 0; i < ny - 2; i++)
		for (int d = 0; d < dim; d++)
			b[(i + 1) * nx * dim + d] = 0.0;

	// x = lx
	for (int i = 0; i < ny - 2; i++)
		for (int d = 0; d < dim; d++)
			b[((i + 2) * nx - 1) * dim + d] = 0.0;

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
void micropp<2>::get_elem_mat(double Ae[2 * 4 * 2 * 4 ], int ex, int ey)
{
	INST_START;

	const int e = glo_elem3D(ex, ey, 0);
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


	for (int i = 0; i < npedim * npedim; i++)
		Ae[i] = 0.0;

	for (int gp = 0; gp < npe; ++gp) {

		double bmat[3][2 * 4], cxb[3][2 * 4] = { 0.0 };
		calc_bmat(gp, bmat);

		for (int i = 0; i < nvoi; ++i) {
			for (int j = 0; j < npedim; ++j) {
				for (int k = 0; k < nvoi; ++k)
					cxb[i][j] += ctan[i][k] * bmat[k][j];
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
	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {
			get_elem_mat(Ae, ex, ey);
			ell_add_struct2D(&A, ex, ey, Ae, dim, nx, ny);
		}
	}
	ell_set_bc_2D(&A, dim, nx, ny);
}


template <>
void micropp<2>::calc_ave_stress(double stress_ave[6])
{
	bool non_linear_flag;
	const double wg = (1 / 4.0) * dx * dy;

	for (int v = 0; v < nvoi; v++)
		stress_ave[v] = 0.0;

	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {

			double stress_aux[3] = { 0.0 };

			for (int gp = 0; gp < npe; ++gp) {

				double stress_gp[6];

				double strain_gp[3];
				get_strain(gp, strain_gp, ex, ey);
				get_stress(gp, strain_gp, &non_linear_flag, stress_gp, ex, ey);
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
void micropp<2>::calc_ave_strain(double strain_ave[6])
{
	const double wg = (1 / 4.0) * dx * dy;

	for (int v = 0; v < nvoi; v++)
		strain_ave[v] = 0.0;

	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {

			double strain_aux[3] = { 0.0 };

			for (int gp = 0; gp < npe; gp++) {
				double strain_gp[6];

				get_strain(gp, strain_gp, ex, ey);
				for (int v = 0; v < nvoi; v++)
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
	bool non_linear_flag;

	const double ivol = 1.0 / (dx * dy);
	const double wg = 0.25 * dx * dy;

	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {

			double strain_aux[3] = { 0.0 };
			double stress_aux[3] = { 0.0 };

			for (int gp = 0; gp < npe; gp++) {

				double stress_gp[3], strain_gp[3];

				get_strain(gp, strain_gp, ex, ey);
				get_stress(gp, strain_gp, &non_linear_flag, stress_gp, ex, ey);
				for (int v = 0; v < nvoi; v++) {
					strain_aux[v] += strain_gp[v] * wg;
					stress_aux[v] += stress_gp[v] * wg;
				}

			}

			const int e = glo_elem3D(ex, ey, 0);
			for (int v = 0; v < nvoi; v++) {
				elem_strain[e * nvoi + v] = strain_aux[v] * ivol;
				elem_stress[e * nvoi + v] = stress_aux[v] * ivol;
			}
		}
	}
}


// Explicit instantiation
template class micropp<2>;
