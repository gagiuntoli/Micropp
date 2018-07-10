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

#include <cmath>
#include <iostream>

#include "instrument.hpp"
#include "micro.hpp"

#define nod_index(i,j,k) ((k)*nx*ny + (j)*nx + (i))

using namespace std;

void micropp_t::set_displ(double *eps)
{

	if (dim == 2) {
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

	} else if (dim == 3) {

		// z = 0
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				const int n = nod_index(i, j, 0);
				const double xcoor = i * dx;
				const double ycoor = j * dy;
				const double zcoor = 0.0;
				const double dux = eps[0] * xcoor + 0.5 * eps[3] * ycoor + 0.5 * eps[4] * zcoor;
				const double duy = 0.5 * eps[3] * xcoor + eps[1] * ycoor + 0.5 * eps[5] * zcoor;
				const double duz = 0.5 * eps[4] * xcoor + 0.5 * eps[5] * ycoor + eps[2] * zcoor;
				u[n * dim] = dux;
				u[n * dim + 1] = duy;
				u[n * dim + 2] = duz;
			}
		}
		// z = lx
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				const int n = nod_index(i, j, nz - 1);
				const double xcoor = i * dx;
				const double ycoor = j * dy;
				const double zcoor = lz;
				const double dux = eps[0] * xcoor + 0.5 * eps[3] * ycoor + 0.5 * eps[4] * zcoor;
				const double duy = 0.5 * eps[3] * xcoor + eps[1] * ycoor + 0.5 * eps[5] * zcoor;
				const double duz = 0.5 * eps[4] * xcoor + 0.5 * eps[5] * ycoor + eps[2] * zcoor;
				u[n * dim] = dux;
				u[n * dim + 1] = duy;
				u[n * dim + 2] = duz;
			}
		}

		// y = 0
		for (int i = 0; i < nx; i++) {
			for (int k = 1; k < nz - 1; k++) {
				const int n = nod_index(i, 0, k);
				const double xcoor = i * dx;
				const double ycoor = 0.0;
				const double zcoor = k * dz;
				const double dux = eps[0] * xcoor + 0.5 * eps[3] * ycoor + 0.5 * eps[4] * zcoor;
				const double duy = 0.5 * eps[3] * xcoor + eps[1] * ycoor + 0.5 * eps[5] * zcoor;
				const double duz = 0.5 * eps[4] * xcoor + 0.5 * eps[5] * ycoor + eps[2] * zcoor;
				u[n * dim] = dux;
				u[n * dim + 1] = duy;
				u[n * dim + 2] = duz;
			}
		}

		// y = ly
		for (int i = 0; i < nx; i++) {
			for (int k = 1; k < nz - 1; k++) {
				const int n = nod_index(i, ny - 1, k);
				const double xcoor = i * dx;
				const double ycoor = ly;
				const double zcoor = k * dz;
				const double dux = eps[0] * xcoor + 0.5 * eps[3] * ycoor + 0.5 * eps[4] * zcoor;
				const double duy = 0.5 * eps[3] * xcoor + eps[1] * ycoor + 0.5 * eps[5] * zcoor;
				const double duz = 0.5 * eps[4] * xcoor + 0.5 * eps[5] * ycoor + eps[2] * zcoor;
				u[n * dim] = dux;
				u[n * dim + 1] = duy;
				u[n * dim + 2] = duz;
			}
		}

		// x=0
		for (int j = 1; j < ny - 1; j++) {
			for (int k = 1; k < nz - 1; k++) {
				const int n = nod_index(0, j, k);
				const double xcoor = 0.0;
				const double ycoor = j * dy;
				const double zcoor = k * dz;
				const double dux = eps[0] * xcoor + 0.5 * eps[3] * ycoor + 0.5 * eps[4] * zcoor;
				const double duy = 0.5 * eps[3] * xcoor + eps[1] * ycoor + 0.5 * eps[5] * zcoor;
				const double duz = 0.5 * eps[4] * xcoor + 0.5 * eps[5] * ycoor + eps[2] * zcoor;
				u[n * dim] = dux;
				u[n * dim + 1] = duy;
				u[n * dim + 2] = duz;
			}
		}

		// x=lx
		for (int j = 1; j < ny - 1; j++) {
			for (int k = 1; k < nz - 1; k++) {
				const int n = nod_index(nx - 1, j, k);
				const double xcoor = lx;
				const double ycoor = j * dy;
				const double zcoor = k * dz;
				const double dux = eps[0] * xcoor + 0.5 * eps[3] * ycoor + 0.5 * eps[4] * zcoor;
				const double duy = 0.5 * eps[3] * xcoor + eps[1] * ycoor + 0.5 * eps[5] * zcoor;
				const double duz = 0.5 * eps[4] * xcoor + 0.5 * eps[5] * ycoor + eps[2] * zcoor;
				u[n * dim] = dux;
				u[n * dim + 1] = duy;
				u[n * dim + 2] = duz;
			}
		}

	}
}

double micropp_t::assembly_rhs(bool *non_linear)
{
	INST_START;

	bool non_linear_one_elem = false;
	*non_linear = false;

	for (int i = 0; i < nn * dim; ++i)
		b[i] = 0.0;

	if (dim == 2) {

		double be[2 * 4];

		for (int ex = 0; ex < nx - 1; ++ex) {
			for (int ey = 0; ey < ny - 1; ++ey) {

				const int n[4] = { ey * nx + ex,
				                   ey * nx + ex + 1,
				                   (ey + 1) * nx + ex + 1,
				                   (ey + 1) * nx + ex };

				const int index[2 * 4] = { n[0] * dim,
				                           n[0] * dim + 1,
				                           n[1] * dim,
				                           n[1] * dim + 1,
				                           n[2] * dim,
				                           n[2] * dim + 1,
				                           n[3] * dim,
				                           n[3] * dim + 1 };

				get_elem_rhs2D(ex, ey, &non_linear_one_elem, be);
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

	} else if (dim == 3) {

		double be[3 * 8];
		int index[3 * 8];

		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {
				for (int ez = 0; ez < nz - 1; ez++) {

					const int n[] = { ez * (nx * ny) + ey * nx + ex,
					                  ez * (nx * ny) + ey * nx + ex + 1,
					                  ez * (nx * ny) + (ey + 1) * nx + ex + 1,
					                  ez * (nx * ny) + (ey + 1) * nx + ex,
					                  n[0] + nx * ny,
					                  n[1] + nx * ny,
					                  n[2] + nx * ny,
					                  n[3] + nx * ny };

					for (int j = 0; j < 8; ++j) {
						for (int d = 0; d < dim; ++d) {
							index[j * dim + d] = n[j] * dim + d;
						}
					}

					get_elem_rhs3D(ex, ey, ez, &non_linear_one_elem, be);
					if (non_linear_one_elem)
						*non_linear = true;

					for (int i = 0; i < npe * dim; ++i)
						b[index[i]] += be[i];

				}
			}
		}

		// boundary conditions
		// z=0
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				const int n = nod_index(i, j, 0);
				for (int d = 0; d < dim; d++)
					b[n * dim + d] = 0.0;
			}
		}
		// z = lx
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				const int n = nod_index(i, j, nz - 1);
				for (int d = 0; d < dim; d++)
					b[n * dim + d] = 0.0;
			}
		}
		// y = 0
		for (int i = 0; i < nx; i++) {
			for (int k = 1; k < nz - 1; k++) {
				const int n = nod_index(i, 0, k);
				for (int d = 0; d < dim; d++)
					b[n * dim + d] = 0.0;
			}
		}
		// y = ly
		for (int i = 0; i < nx; i++) {
			for (int k = 1; k < nz - 1; k++) {
				const int n = nod_index(i, ny - 1, k);
				for (int d = 0; d < dim; d++)
					b[n * dim + d] = 0.0;
			}
		}
		// x=0
		for (int j = 1; j < ny - 1; j++) {
			for (int k = 1; k < nz - 1; k++) {
				const int n = nod_index(0, j, k);
				for (int d = 0; d < dim; d++)
					b[n * dim + d] = 0.0;
			}
		}
		// x=lx
		for (int j = 1; j < ny - 1; j++) {
			for (int k = 1; k < nz - 1; k++) {
				const int n = nod_index(nx - 1, j, k);
				for (int d = 0; d < dim; d++)
					b[n * dim + d] = 0.0;
			}
		}
	}

	for (int i = 0; i < nn * dim; ++i)
		b[i] = -b[i];

	double norm = 0.0;
	for (int i = 0; i < nn * dim; ++i)
		norm += b[i] * b[i];
	norm = sqrt(norm);

	return norm;
}

void micropp_t::assembly_mat()
{
	INST_START;

	int index[8];

	ell_set_zero_mat(&A);

	if (dim == 2) {

		double Ae[2 * 4 * 2 * 4];
		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {
				get_elem_mat2D(ex, ey, Ae);
				ell_add_struct2D(&A, ex, ey, Ae, dim, nx, ny);
			}
		}
		ell_set_bc_2D(&A, dim, nx, ny);

	} else if (dim == 3) {

		double Ae[3 * 8 * 3 * 8];
		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {
				for (int ez = 0; ez < nz - 1; ez++) {
					get_elem_mat3D(ex, ey, ez, Ae);
					ell_add_struct3D(&A, ex, ey, ez, Ae, dim, nx, ny, nz);
				}
			}
		}
		ell_set_bc_3D(&A, dim, nx, ny, nz);
	}
}

void micropp_t::get_elem_mat2D(int ex, int ey, double (&Ae)[2 * 4 * 2 * 4])
{
	INST_START;

	double nu, E;
	double ctan[3][3];
	bool plasticity;

	int e = glo_elem3D(ex, ey, 0);

	const material_t material = get_material(e);

	E = material.E;
	nu = material.nu;
	plasticity = material.plasticity;

	ctan[0][0] = (1 - nu);
	ctan[0][1] = nu;
	ctan[0][2] = 0;
	ctan[1][0] = nu;
	ctan[1][1] = (1 - nu);
	ctan[1][2] = 0;
	ctan[2][0] = 0;
	ctan[2][1] = 0;
	ctan[2][2] = (1 - 2 * nu) / 2;

	for (int i = 0; i < nvoi; i++)
		for (int j = 0; j < nvoi; j++)
			ctan[i][j] *= E / ((1 + nu) * (1 - 2 * nu));

	double xg[4][2] = {
		{-0.577350269189626, -0.577350269189626},
		{+0.577350269189626, -0.577350269189626},
		{+0.577350269189626, +0.577350269189626},
		{-0.577350269189626, +0.577350269189626}
	};

	double dsh[4][2], b_mat[3][8], cxb[3][8];

	for (int i = 0; i < npe * dim * npe * dim; i++)
		Ae[i] = 0.0;

	for (int gp = 0; gp < 4; gp++) {

		dsh[0][0] = -(1 - xg[gp][1]) / 4 * 2 / dx;
		dsh[0][1] = -(1 - xg[gp][0]) / 4 * 2 / dy;
		dsh[1][0] = +(1 - xg[gp][1]) / 4 * 2 / dx;
		dsh[1][1] = -(1 + xg[gp][0]) / 4 * 2 / dy;
		dsh[2][0] = +(1 + xg[gp][1]) / 4 * 2 / dx;
		dsh[2][1] = +(1 + xg[gp][0]) / 4 * 2 / dy;
		dsh[3][0] = -(1 + xg[gp][1]) / 4 * 2 / dx;
		dsh[3][1] = +(1 - xg[gp][0]) / 4 * 2 / dy;

		for (int i = 0; i < 4; i++) {
			b_mat[0][i * dim] = dsh[i][0];
			b_mat[0][i * dim + 1] = 0;
			b_mat[1][i * dim] = 0;
			b_mat[1][i * dim + 1] = dsh[i][1];
			b_mat[2][i * dim] = dsh[i][1];
			b_mat[2][i * dim + 1] = dsh[i][0];
		}

		for (int i = 0; i < nvoi; i++) {
			for (int j = 0; j < npe * dim; j++) {
				cxb[i][j] = 0.0;
				for (int k = 0; k < nvoi; k++)
					cxb[i][j] += ctan[i][k] * b_mat[k][j];
			}
		}

		double wg = 0.25 * dx * dy;
		for (int i = 0; i < npe * dim; i++)
			for (int j = 0; j < npe * dim; j++)
				for (int m = 0; m < nvoi; m++)
					Ae[i * npe * dim + j] += b_mat[m][i] * cxb[m][j] * wg;

	}			// gp loop
}

void micropp_t::get_elem_mat3D(int ex, int ey, int ez, double (&Ae)[3 * 8 * 3 * 8])
{
	INST_START;
	const int e = glo_elem3D(ex, ey, ez);
	const material_t material = get_material(e);
	double ctan[6][6];
	const int npedim = npe * dim;
	const double wg = (1 / 8.0) * dx * dy * dz;

	/*
	  C = lambda * (1x1) + 2 mu I
	  last 3 components are without the 2 because we use eps = {e11 e22 e33 2*e12 2*e13 2*e23}
	*/

	for (int i = 0; i < npe * dim * npe * dim; i++)
		Ae[i] = 0.0;

	for (int gp = 0; gp < 8; gp++) {

		if (material.plasticity) {
			get_ctan_plast_pert(ex, ey, ez, gp, ctan);   // ctan_pert
			//get_ctan_plast_sec(ex, ey, ez, gp, ctan);    // ctan_secant
			//get_ctan_plast_exact(ex, ey, ez, gp, ctan);  // ctan_exact
		} else {

			for (int i = 0; i < 6; ++i)
				for (int j = 0; j < 6; ++j)
					ctan[i][j] = 0.0;

			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					ctan[i][j] += material.lambda;

			for (int i = 0; i < 3; ++i)
				ctan[i][i] += 2 * material.mu;

			for (int i = 3; i < 6; ++i)
				ctan[i][i] += material.mu;

		}

		double bmat[6][3 * 8], cxb[6][3 * 8] = { 0.0 };
		calc_bmat_3D(gp, bmat);

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
	}			// gp loop
}

void micropp_t::get_ctan_plast_sec(int ex, int ey, int ez, int gp, double ctan[6][6])
{
	INST_START;
	bool non_linear;
	double stress_pert[6], strain_pert[6], strain_0[6], d_strain = 1.0e-8;
	get_strain3D(ex, ey, ez, gp, strain_0);

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++)
			strain_pert[j] = 0.0;

		if (fabs(strain_0[i]) > 1.0e-7)
			strain_pert[i] = strain_0[i];
		else
			strain_pert[i] = d_strain;

		get_stress3D(ex, ey, ez, gp, strain_pert, &non_linear, stress_pert);

		for (int j = 0; j < 6; j++)
			ctan[j][i] = stress_pert[j] / strain_pert[i];
	}

}

void micropp_t::get_ctan_plast_exact(int ex, int ey, int ez, int gp, double ctan[6][6])
{
	INST_START;
	double strain[6];
	int e = glo_elem3D(ex, ey, ez);
	get_strain3D(ex, ey, ez, gp, strain);

	const material_t material = get_material(e);

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			ctan[i][j] = 0.0;

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ctan[i][j] += material.lambda;

	for (int i = 0; i < 3; i++)
		ctan[i][i] += 2 * material.mu;

	for (int i = 3; i < 6; i++)
		ctan[i][i] += material.mu;

	//  double theta_1 = 1 - 2*material.mu*dl / sig_dev_trial_norm;
	//  double theta_2 = 1 / (1 + material.Ka) - (1 - theta_1);
	//  for (int i=0; i<6; i++)
	//    for (int j=0; j<6; j++)
	//      ctan[i][j] -= 2 * material.mu * theta_2 * normal[i] * normal[j];
	//
	//  for (int i=0; i<3; i++)
	//    for (int j=0; j<3; j++)
	//      ctan[i][j] -= 2 * material.mu * theta_1 * (1.0/3);
	//
	//  for (int i=0; i<3; i++)
	//    ctan[i][i] += 2 * material.mu * theta_1;
}

void micropp_t::get_ctan_plast_pert(int ex, int ey, int ez, int gp, double ctan[6][6])
{
	INST_START;
	bool non_linear;
	double stress_0[6], stress_pert[6], strain_0[6], strain_pert[6], deps = 1.0e-8;
	get_strain3D(ex, ey, ez, gp, strain_0);
	get_stress3D(ex, ey, ez, gp, strain_0, &non_linear, stress_0);

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++)
			strain_pert[j] = strain_0[j];

		strain_pert[i] += deps;
		get_stress3D(ex, ey, ez, gp, strain_pert, &non_linear, stress_pert);

		for (int j = 0; j < 6; j++)
			ctan[j][i] = (stress_pert[j] - stress_0[j]) / deps;
	}
}

void micropp_t::calc_bmat_3D(int gp, double bmat[6][3 * 8])
{
	const double dsh[8][3] = {
		{ -(1 - xg[gp][1]) * (1 - xg[gp][2]) * 0.25 / dx,
		  -(1 - xg[gp][0]) * (1 - xg[gp][2]) * 0.25 / dy,
		  -(1 - xg[gp][0]) * (1 - xg[gp][1]) * 0.25 / dz },
		{ +(1 - xg[gp][1]) * (1 - xg[gp][2]) * 0.25 / dx,
		  -(1 + xg[gp][0]) * (1 - xg[gp][2]) * 0.25 / dy,
		  -(1 + xg[gp][0]) * (1 - xg[gp][1]) * 0.25 / dz },
		{ +(1 + xg[gp][1]) * (1 - xg[gp][2]) * 0.25 / dx,
		  +(1 + xg[gp][0]) * (1 - xg[gp][2]) * 0.25 / dy,
		  -(1 + xg[gp][0]) * (1 + xg[gp][1]) * 0.25 / dz },
		{ -(1 + xg[gp][1]) * (1 - xg[gp][2]) * 0.25 / dx,
		  +(1 - xg[gp][0]) * (1 - xg[gp][2]) * 0.25 / dy,
		  -(1 - xg[gp][0]) * (1 + xg[gp][1]) * 0.25 / dz },
		{ -(1 - xg[gp][1]) * (1 + xg[gp][2]) * 0.25 / dx,
		  -(1 - xg[gp][0]) * (1 + xg[gp][2]) * 0.25 / dy,
		  +(1 - xg[gp][0]) * (1 - xg[gp][1]) * 0.25 / dz },
		{ +(1 - xg[gp][1]) * (1 + xg[gp][2]) * 0.25 / dx,
		  -(1 + xg[gp][0]) * (1 + xg[gp][2]) * 0.25 / dy,
		  +(1 + xg[gp][0]) * (1 - xg[gp][1]) * 0.25 / dz },
		{ +(1 + xg[gp][1]) * (1 + xg[gp][2]) * 0.25 / dx,
		  +(1 + xg[gp][0]) * (1 + xg[gp][2]) * 0.25 / dy,
		  +(1 + xg[gp][0]) * (1 + xg[gp][1]) * 0.25 / dz },
		{ -(1 + xg[gp][1]) * (1 + xg[gp][2]) * 0.25 / dx,
		  +(1 - xg[gp][0]) * (1 + xg[gp][2]) * 0.25 / dy,
		  +(1 - xg[gp][0]) * (1 + xg[gp][1]) * 0.25 / dz } };

	for (int i = 0; i < 8; i++) {
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

void micropp_t::get_elem_rhs2D(int ex, int ey, bool *non_linear, double (&be)[2 * 4])
{
	double bmat[3][2 * 4], cxb[3][8], stress_gp[6];
	const double wg = 0.25 * dx * dy;

	for (int i = 0; i < 2 * 4; i++)
		be[i] = 0.0;

	for (int gp = 0; gp < 4; gp++) {

		const double dsh[4][2] = { { -(1 - xg[gp][1]) * 0.5 / dx,
		                             -(1 - xg[gp][0]) * 0.5 / dy },
		                           { +(1 - xg[gp][1]) * 0.5 / dx,
		                             -(1 + xg[gp][0]) * 0.5 / dy },
		                           { +(1 + xg[gp][1]) * 0.5 / dx,
		                             +(1 + xg[gp][0]) * 0.5 / dy },
		                           { -(1 + xg[gp][1]) * 0.5 / dx,
		                             +(1 - xg[gp][0]) * 0.5 / dy } };

		for (int i = 0; i < 4; i++) {
			bmat[0][i * dim] = dsh[i][0];
			bmat[0][i * dim + 1] = 0;
			bmat[1][i * dim] = 0;
			bmat[1][i * dim + 1] = dsh[i][1];
			bmat[2][i * dim] = dsh[i][1];
			bmat[2][i * dim + 1] = dsh[i][0];
		}

		double strain_gp[3];
		get_strain2D(ex, ey, gp, strain_gp);
		get_stress2D(ex, ey, gp, strain_gp, non_linear, stress_gp);

		for (int i = 0; i < npe * dim; i++) {
			for (int j = 0; j < nvoi; j++) {
				be[i] += bmat[j][i] * stress_gp[j] * wg;
			}
		}

	}			// gp loop
}

void micropp_t::get_elem_rhs3D(int ex, int ey, int ez, bool * non_linear, double (&be)[3 * 8])
{
	double bmat[6][3 * 8], cxb[6][3 * 8], stress_gp[6];

	for (int i = 0; i < 3 * 8; i++)
		be[i] = 0.0;

	for (int gp = 0; gp < 8; gp++) {

		calc_bmat_3D(gp, bmat);

		double strain_gp[6];
		get_strain3D(ex, ey, ez, gp, strain_gp);
		get_stress3D(ex, ey, ez, gp, strain_gp, non_linear, stress_gp);

		double wg = (1 / 8.0) * dx * dy * dz;
		for (int i = 0; i < npe * dim; i++)
			for (int j = 0; j < nvoi; j++)
				be[i] += bmat[j][i] * stress_gp[j] * wg;

	}			// gp loop
}

void micropp_t::calc_ave_stress(double stress_ave[6])
{
	bool non_linear_flag;

	for (int v = 0; v < nvoi; v++)
		stress_ave[v] = 0.0;

	if (dim == 2) {

		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {

				double stress_aux[3];
				for (int i = 0; i < nvoi; i++)
					stress_aux[i] = 0.0;

				for (int gp = 0; gp < 4; gp++) {

					double stress_gp[6];
					double wg = 0.25 * dx * dy;

					double strain_gp[3];
					get_strain2D(ex, ey, gp, strain_gp);
					get_stress2D(ex, ey, gp, strain_gp, &non_linear_flag, stress_gp);
					for (int v = 0; v < nvoi; v++)
						stress_aux[v] += stress_gp[v] * wg;

				}
				for (int v = 0; v < nvoi; v++)
					stress_ave[v] += stress_aux[v];
			}
		}

	} else if (dim == 3) {

		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {
				for (int ez = 0; ez < nz - 1; ez++) {

					double stress_aux[6];
					for (int i = 0; i < nvoi; i++)
						stress_aux[i] = 0.0;

					for (int gp = 0; gp < 8; gp++) {

						double stress_gp[6];
						double wg = (1 / 8.0) * dx * dy * dz;

						double strain_gp[6];
						get_strain3D(ex, ey, ez, gp, strain_gp);
						get_stress3D(ex, ey, ez, gp, strain_gp, &non_linear_flag, stress_gp);
						for (int v = 0; v < nvoi; v++)
							stress_aux[v] += stress_gp[v] * wg;

					}
					for (int v = 0; v < nvoi; v++)
						stress_ave[v] += stress_aux[v];
				}
			}
		}
	}

	for (int v = 0; v < nvoi; v++)
		stress_ave[v] /= (lx * ly);
}

void micropp_t::calc_ave_strain(double strain_ave[6])
{
	for (int v = 0; v < nvoi; v++)
		strain_ave[v] = 0.0;

	if (dim == 2) {

		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {

				double strain_aux[3];
				for (int i = 0; i < nvoi; i++)
					strain_aux[i] = 0.0;

				for (int gp = 0; gp < 4; gp++) {
					double strain_gp[6];
					double wg = 0.25 * dx * dy;
					get_strain2D(ex, ey, gp, strain_gp);
					for (int v = 0; v < nvoi; v++)
						strain_aux[v] += strain_gp[v] * wg;
				}

				for (int v = 0; v < nvoi; v++)
					strain_ave[v] += strain_aux[v];
			}
		}

	} else if (dim == 3) {

		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {
				for (int ez = 0; ez < nz - 1; ez++) {

					double strain_aux[6];
					for (int i = 0; i < nvoi; i++)
						strain_aux[i] = 0.0;

					for (int gp = 0; gp < 8; gp++) {
						double strain_gp[6];
						double wg = (1 / 8.0) * dx * dy * dz;
						get_strain3D(ex, ey, ez, gp, strain_gp);
						for (int v = 0; v < nvoi; v++)
							strain_aux[v] += strain_gp[v] * wg;
					}

					for (int v = 0; v < nvoi; v++)
						strain_ave[v] += strain_aux[v];
				}
			}
		}
	}

	for (int v = 0; v < nvoi; v++)
		strain_ave[v] /= (lx * ly);
}

void micropp_t::calc_fields()
{
	bool non_linear_flag;

	if (dim == 2) {

		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {

				double strain_aux[3];
				double stress_aux[3];
				for (int i = 0; i < nvoi; i++) {
					strain_aux[i] = 0.0;
					stress_aux[i] = 0.0;
				}

				for (int gp = 0; gp < 4; gp++) {

					double stress_gp[3], strain_gp[3];
					double wg = 0.25 * dx * dy;

					get_strain2D(ex, ey, gp, strain_gp);
					get_stress2D(ex, ey, gp, strain_gp, &non_linear_flag, stress_gp);
					for (int v = 0; v < nvoi; v++) {
						strain_aux[v] += strain_gp[v] * wg;
						stress_aux[v] += stress_gp[v] * wg;
					}

				}
				double vol = dx * dy;
				int e = glo_elem3D(ex, ey, 0);
				for (int v = 0; v < nvoi; v++) {
					elem_strain[e * nvoi + v] = strain_aux[v] / vol;
					elem_stress[e * nvoi + v] = stress_aux[v] / vol;
				}
			}
		}

	} else if (dim == 3) {

		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {
				for (int ez = 0; ez < nz - 1; ez++) {

					double strain_aux[6];
					double stress_aux[6];
					for (int i = 0; i < nvoi; i++) {
						strain_aux[i] = 0.0;
						stress_aux[i] = 0.0;
					}

					for (int gp = 0; gp < 8; gp++) {

						double stress_gp[6], strain_gp[6];
						double wg = (1 / 8.0) * dx * dy * dz;

						get_strain3D(ex, ey, ez, gp, strain_gp);
						get_stress3D(ex, ey, ez, gp, strain_gp, &non_linear_flag, stress_gp);
						for (int v = 0; v < nvoi; v++) {
							strain_aux[v] += strain_gp[v] * wg;
							stress_aux[v] += stress_gp[v] * wg;
						}

					}

					double vol = dx * dy * dz;
					int e = glo_elem3D(ex, ey, ez);
					for (int v = 0; v < nvoi; v++) {
						elem_strain[e * nvoi + v] = strain_aux[v] / vol;
						elem_stress[e * nvoi + v] = stress_aux[v] / vol;
					}
				}
			}
		}

	}
}

void micropp_t::get_stress2D(int ex, int ey, int gp, double strain_gp[3], bool * non_linear, double *stress_gp)
{
	*non_linear = false;

	double nu, E;
	double ctan[3][3];
	bool plasticity;

	int e = glo_elem3D(ex, ey, 0);

	const material_t material = get_material(e);

	E = material.E;
	nu = material.nu;
	plasticity = material.plasticity;

	ctan[0][0] = (1 - nu);
	ctan[0][1] = nu;
	ctan[0][2] = 0;
	ctan[1][0] = nu;
	ctan[1][1] = (1 - nu);
	ctan[1][2] = 0;
	ctan[2][0] = 0;
	ctan[2][1] = 0;
	ctan[2][2] = (1 - 2 * nu) / 2;
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

void micropp_t::get_dev_tensor(double tensor[6], double tensor_dev[6])
{
	for (int i = 0; i < 6; i++)
		tensor_dev[i] = tensor[i];
	for (int i = 0; i < 3; i++)
		tensor_dev[i] -= (1 / 3.0) * (tensor[0] + tensor[1] + tensor[2]);

}

void micropp_t::get_stress3D(int ex, int ey, int ez, int gp, double eps[6],
                             bool * non_linear, double *stress_gp)
{
	*non_linear = false;

	double ctan[6][6];
	int e = glo_elem3D(ex, ey, ez);
	const material_t material = get_material(e);

	if (material.plasticity == true) {
		double alpha_old, alpha_new, eps_p_old[6], eps_p_new[6];
		for (int i = 0; i < 6; ++i)
			eps_p_old[i] = vars_old[intvar_ix(e, gp, i)];
		alpha_old = vars_old[intvar_ix(e, gp, 6)];

		plastic_step(&material, eps, eps_p_old, alpha_old, eps_p_new,
		             &alpha_new, non_linear, stress_gp);

		for (int i = 0; i < 6; ++i)
			vars_new[intvar_ix(e, gp, i)] = eps_p_new[i];
		vars_new[intvar_ix(e, gp, 6)] = alpha_new;

	} else {

		for (int i = 0; i < 3; ++i) {
			stress_gp[i] = material.lambda * (eps[0] + eps[1] + eps[2]);
			stress_gp[i] += 2 * material.mu * eps[i];
		}
		for (int i = 3; i < 6; ++i)
			stress_gp[i] = material.mu * eps[i];
	}

}

void micropp_t::plastic_step(const material_t *material, double eps[6],
                             double eps_p_old[6], double alpha_old,
                             double eps_p_new[6], double *alpha_new,
                             bool *non_linear, double stress[6])
{
	double eps_dev[6];
	double eps_p_dev_1[6];
	double normal[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double dl = 0.0;
	double sig_dev_trial[6], sig_dev_trial_norm;

	get_dev_tensor(eps_p_old, eps_p_dev_1);
	get_dev_tensor(eps, eps_dev);

	for (int i = 0; i < 3; ++i)
		sig_dev_trial[i] = 2 * material->mu * (eps_dev[i] - eps_p_dev_1[i]);

	for (int i = 3; i < 6; ++i)
		sig_dev_trial[i] = material->mu * (eps_dev[i] - eps_p_dev_1[i]);

	sig_dev_trial_norm = sqrt(sig_dev_trial[0] * sig_dev_trial[0] +
	                          sig_dev_trial[1] * sig_dev_trial[1] +
	                          sig_dev_trial[2] * sig_dev_trial[2] +
	                          2 * sig_dev_trial[3] * sig_dev_trial[3] +
	                          2 * sig_dev_trial[4] * sig_dev_trial[4] +
	                          2 * sig_dev_trial[5] * sig_dev_trial[5]);

	double f_trial = sig_dev_trial_norm - sqrt(2.0 / 3) * (material->Sy + material->Ka * alpha_old);

	if (f_trial > 0) {
		*non_linear = true;

		for (int i = 0; i < 6; ++i)
			normal[i] = sig_dev_trial[i] / sig_dev_trial_norm;

		dl = f_trial / (2 * material->mu * (1.0 + (0.0 * material->Ka) / (3 * material->mu)));

		for (int i = 0; i < 6; ++i)
			eps_p_new[i] = eps_p_old[i] + dl * normal[i];
		*alpha_new = alpha_old + sqrt(2.0 / 3) * dl;
	} else {
		for (int i = 0; i < 6; ++i)
			eps_p_new[i] = eps_p_old[i];
		*alpha_new = alpha_old;
	}

	//sig_2 = s_trial + K * tr(eps) * 1 - 2*mu*dl*normal;
	for (int i = 0; i < 6; ++i)
		stress[i] = sig_dev_trial[i];

	for (int i = 0; i < 3; ++i)
		stress[i] += material->k * (eps[0] + eps[1] + eps[2]);

	for (int i = 0; i < 6; ++i)
		stress[i] -= 2 * material->mu * dl * normal[i];
}

material_t micropp_t::get_material(const int e)
{
	int mat_num;
	if (micro_type == 0)
		if (elem_type[e] == 0)
			mat_num = 0;
		else
			mat_num = 1;

	else if (micro_type == 1)
		if (elem_type[e] == 0)
			mat_num = 0;
		else
			mat_num = 1;

	return material_list[mat_num];
}

void micropp_t::get_strain2D(int ex, int ey, int gp, double *strain_gp)
{
	double elem_disp[2 * 4];
	getElemDisp(ex, ey, elem_disp);

	double xg[4][2] = {
		{-0.577350269189626, -0.577350269189626},
		{+0.577350269189626, -0.577350269189626},
		{+0.577350269189626, +0.577350269189626},
		{-0.577350269189626, +0.577350269189626}
	};

	double dsh[4][2], bmat[3][8], cxb[3][8];

	dsh[0][0] = -(1 - xg[gp][1]) / 4 * 2 / dx;
	dsh[0][1] = -(1 - xg[gp][0]) / 4 * 2 / dy;
	dsh[1][0] = +(1 - xg[gp][1]) / 4 * 2 / dx;
	dsh[1][1] = -(1 + xg[gp][0]) / 4 * 2 / dy;
	dsh[2][0] = +(1 + xg[gp][1]) / 4 * 2 / dx;
	dsh[2][1] = +(1 + xg[gp][0]) / 4 * 2 / dy;
	dsh[3][0] = -(1 + xg[gp][1]) / 4 * 2 / dx;
	dsh[3][1] = +(1 - xg[gp][0]) / 4 * 2 / dy;

	for (int i = 0; i < 4; i++) {
		bmat[0][i * dim] = dsh[i][0];
		bmat[0][i * dim + 1] = 0;
		bmat[1][i * dim] = 0;
		bmat[1][i * dim + 1] = dsh[i][1];
		bmat[2][i * dim] = dsh[i][1];
		bmat[2][i * dim + 1] = dsh[i][0];
	}
	for (int v = 0; v < nvoi; v++) {
		strain_gp[v] = 0.0;
		for (int i = 0; i < npe * dim; i++) {
			strain_gp[v] += bmat[v][i] * elem_disp[i];
		}
	}

}

void micropp_t::get_strain3D(int ex, int ey, int ez, int gp, double *strain_gp)
{

	double elem_disp[3 * 8];
	getElemDisp(ex, ey, ez, elem_disp);

	double bmat[6][3 * 8];
	calc_bmat_3D(gp, bmat);

	for (int v = 0; v < nvoi; ++v) {
		strain_gp[v] = 0.0;
		for (int i = 0; i < npe * dim; i++)
			strain_gp[v] += bmat[v][i] * elem_disp[i];
	}

}

void micropp_t::getElemDisp(int ex, int ey, double *elem_disp)
{
	int n0 = ey * nx + ex;
	int n1 = ey * nx + ex + 1;
	int n2 = (ey + 1) * nx + ex + 1;
	int n3 = (ey + 1) * nx + ex;

	for (int d = 0; d < dim; d++) {
		elem_disp[0 * dim + d] = u[n0 * dim + d];
		elem_disp[1 * dim + d] = u[n1 * dim + d];
		elem_disp[2 * dim + d] = u[n2 * dim + d];
		elem_disp[3 * dim + d] = u[n3 * dim + d];
	}
}

void micropp_t::getElemDisp(int ex, int ey, int ez, double *elem_disp)
{
	int n0 = ez * (nx * ny) + ey * nx + ex;
	int n1 = ez * (nx * ny) + ey * nx + ex + 1;
	int n2 = ez * (nx * ny) + (ey + 1) * nx + ex + 1;
	int n3 = ez * (nx * ny) + (ey + 1) * nx + ex;
	int n4 = n0 + nx * ny;
	int n5 = n1 + nx * ny;
	int n6 = n2 + nx * ny;
	int n7 = n3 + nx * ny;

	for (int d = 0; d < dim; d++) {
		elem_disp[0 * dim + d] = u[n0 * dim + d];
		elem_disp[1 * dim + d] = u[n1 * dim + d];
		elem_disp[2 * dim + d] = u[n2 * dim + d];
		elem_disp[3 * dim + d] = u[n3 * dim + d];
		elem_disp[4 * dim + d] = u[n4 * dim + d];
		elem_disp[5 * dim + d] = u[n5 * dim + d];
		elem_disp[6 * dim + d] = u[n6 * dim + d];
		elem_disp[7 * dim + d] = u[n7 * dim + d];
	}
}
