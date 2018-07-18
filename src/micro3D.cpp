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
int micropp<3>::get_elem_type(int ex, int ey, int ez)
{
	assert(micro_type == 0 || micro_type == 1);

	if (micro_type == 0) {
		// esfera en matriz
		const double x1 = ex * dx + dx / 2;
		const double y1 = ey * dy + dy / 2;
		const double z1 = ez * dz + dz / 2;
		const double x2 = lx / 2;
		const double y2 = ly / 2;
		const double z2 = lz / 2;
		const double rad = micro_params[3];
		return ((x2 - x1) * (x2 - x1) +
		        (y2 - y1) * (y2 - y1) +
		        (z2 - z1) * (z2 - z1) < width * width);

	} else if (micro_type == 1) {
		const double y = ey * dy + dy / 2;
		return (y < width);
	}

	cerr << "Invalid micro_type = " << micro_type << endl;
	return -1;
}


template<>
micropp<3>::micropp(const int _ngp, const int size[3], const int _micro_type,
                    const double *_micro_params, const int *_mat_types,
                    const double *_params):
	ngp(_ngp),
	nx(size[0]), ny(size[1]), nz(size[2]),
	nn(nx * ny * nz),
	nex(nx - 1), ney(ny - 1), nez(nz - 1),
	nelem(nex * ney * nez),
	lx(_micro_params[0]), ly(_micro_params[1]),	lz(_micro_params[2]),
	dx(lx / nex), dy(ly / ney),	dz(lz / nez),

	width(_micro_params[3]), inv_tol(micro_params[4]),
	micro_type(_micro_type), num_int_vars(nelem * 8 * NUM_VAR_GP)
{
	assert(dim == 3);

	initialize(_micro_params, _mat_types, _params);

	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {
			for (int ez = 0; ez < nz - 1; ez++) {
				int e = glo_elem3D(ex, ey, ez);
				elem_type[e] = get_elem_type(ex, ey, ez);
			}
		}
	}
	ell_init_3D(&A, dim, nx, ny, nz);

	calc_ctan_lin();
}


template <>
void micropp<3>::set_displ(double *eps)
{
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


template <>
template <>
void micropp<3>::get_stress(int gp, double eps[6], bool * non_linear,
                            double *stress_gp, int ex, int ey, int ez) const
{
	*non_linear = false;

	double ctan[6][6];
	const int e = glo_elem3D(ex, ey, ez);
	const material_t material = get_material(e);
	const double mu = material.mu;

	if (material.plasticity == true) {
		double *eps_p_new = NULL, * alpha_new = NULL;
		double * const eps_p_old = &vars_old[intvar_ix(e, gp, 0)];
		const double alpha_old = vars_old[intvar_ix(e, gp, 6)];

		if (vars_new) {
			eps_p_new = &vars_new[intvar_ix(e, gp, 0)];
			alpha_new = &vars_new[intvar_ix(e, gp, 6)];
		}

		plastic_step(&material, eps, eps_p_old, alpha_old, eps_p_new,
		             alpha_new, non_linear, stress_gp);

	} else {

		for (int i = 0; i < 3; ++i)
			stress_gp[i] = material.lambda * (eps[0] + eps[1] + eps[2]) \
						   + 2 * mu * eps[i];

		for (int i = 3; i < 6; ++i)
			stress_gp[i] = mu * eps[i];
	}

}


template <>
template <>
void micropp<3>::calc_bmat(int gp, double bmat[6][3 * 8]) const
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


template <>
template <>
void micropp<3>::getElemDisp(double *elem_disp, int ex, int ey, int ez) const
{
	const int nxny = nx * ny;
	const int n[] = { ez * nxny + ey * nx + ex,
	                  ez * nxny + ey * nx + ex + 1,
	                  ez * nxny + (ey + 1) * nx + ex + 1,
	                  ez * (nx * ny) + (ey + 1) * nx + ex,
	                  n[0] + nxny,
	                  n[1] + nxny,
	                  n[2] + nxny,
	                  n[3] + nxny };

	for (int i = 0 ; i < npe; ++i)
		for (int d = 0; d < dim; ++d)
			elem_disp[i * dim + d] = u[n[i] * dim + d];
}


template <>
template <>
void micropp<3>::get_strain(int gp, double *strain_gp,
                            int ex, int ey, int ez) const
{

	double elem_disp[3 * 8];
	getElemDisp(elem_disp, ex, ey, ez);

	double bmat[6][3 * 8];
	calc_bmat(gp, bmat);

	for (int v = 0; v < nvoi; ++v) {
		strain_gp[v] = 0.0;
		for (int i = 0; i < npe * dim; i++)
			strain_gp[v] += bmat[v][i] * elem_disp[i];
	}
}


template <>
template <>
void micropp<3>::get_elem_rhs(bool *non_linear,
                             double *be, int ex, int ey, int ez) const
{
	double bmat[6][3 * 8], cxb[6][3 * 8], stress_gp[6];
	const double wg = (1 / 8.0) * dx * dy * dz;

	for (int i = 0; i < 3 * 8; i++)
		be[i] = 0.0;

	for (int gp = 0; gp < 8; gp++) {

		calc_bmat(gp, bmat);

		double strain_gp[6];
		get_strain(gp, strain_gp, ex, ey, ez);
		get_stress(gp, strain_gp, non_linear, stress_gp, ex, ey, ez);

		for (int i = 0; i < npe * dim; i++)
			for (int j = 0; j < nvoi; j++)
				be[i] += bmat[j][i] * stress_gp[j] * wg;

	}			// gp loop
}


template<>
double micropp<3>::assembly_rhs(bool *non_linear)
{
	INST_START;

	bool non_linear_one_elem = false;
	*non_linear = false;

	for (int i = 0; i < nn * dim; ++i)
		b[i] = 0.0;

	double be[3 * 8];
	int index[3 * 8];

	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {
			for (int ez = 0; ez < nz - 1; ez++) {

				const int n[npe] = { ez * (nx * ny) + ey * nx + ex,
				                     ez * (nx * ny) + ey * nx + ex + 1,
				                     ez * (nx * ny) + (ey + 1) * nx + ex + 1,
				                     ez * (nx * ny) + (ey + 1) * nx + ex,
				                     n[0] + nx * ny,
				                     n[1] + nx * ny,
				                     n[2] + nx * ny,
				                     n[3] + nx * ny };

				for (int j = 0; j < npe; ++j) {
					for (int d = 0; d < dim; ++d) {
						index[j * dim + d] = n[j] * dim + d;
					}
				}

				get_elem_rhs(&non_linear_one_elem, be, ex, ey, ez);
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

	// Common
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
void micropp<3>::get_elem_mat(double Ae[3 * 8 * 3 * 8], int ex, int ey, int ez)
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

	for (int i = 0; i < npedim * npedim; i++)
		Ae[i] = 0.0;

	for (int gp = 0; gp < npe; ++gp) {

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
void micropp<3>::assembly_mat()
{
	INST_START;

	ell_set_zero_mat(&A);

	double Ae[3 * 8 * 3 * 8];
	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {
			for (int ez = 0; ez < nz - 1; ez++) {
				get_elem_mat(Ae, ex, ey, ez);
				ell_add_struct3D(&A, ex, ey, ez, Ae, dim, nx, ny, nz);
			}
		}
	}
	ell_set_bc_3D(&A, dim, nx, ny, nz);
}


template <>
void micropp<3>::calc_ave_stress(double stress_ave[6])
{
	bool non_linear_flag;
	const double wg = (1 / 8.0) * dx * dy * dz;

	for (int v = 0; v < nvoi; v++)
		stress_ave[v] = 0.0;

	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {
			for (int ez = 0; ez < nz - 1; ez++) {

				double stress_aux[6] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[6];

					double strain_gp[6];
					get_strain(gp, strain_gp, ex, ey, ez);
					get_stress(gp, strain_gp, &non_linear_flag, stress_gp,
					           ex, ey, ez);
					for (int v = 0; v < nvoi; ++v)
						stress_aux[v] += stress_gp[v] * wg;

				}
				for (int v = 0; v < nvoi; ++v)
					stress_ave[v] += stress_aux[v];
			}
		}
	}

	for (int v = 0; v < nvoi; ++v)
		stress_ave[v] /= (lx * ly);
}


template <>
void micropp<3>::calc_ave_strain(double strain_ave[6])
{
	const double wg = (1 / 8.0) * dx * dy * dz;

	for (int v = 0; v < nvoi; v++)
		strain_ave[v] = 0.0;

	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {
			for (int ez = 0; ez < nz - 1; ez++) {

				double strain_aux[6] = { 0.0 };

				for (int gp = 0; gp < npe; gp++) {
					double strain_gp[6];

					get_strain(gp, strain_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; v++)
						strain_aux[v] += strain_gp[v] * wg;
				}

				for (int v = 0; v < nvoi; v++)
					strain_ave[v] += strain_aux[v];
			}
		}
	}


	for (int v = 0; v < nvoi; v++)
		strain_ave[v] /= (lx * ly);
}


template<>
void micropp<3>::calc_fields()
{
	bool non_linear_flag;

	const double ivol = 1.0 / (dx * dy * dz);
	const double wg = (1 / 8.0) * dx * dy * dz;

	for (int ex = 0; ex < nx - 1; ex++) {
		for (int ey = 0; ey < ny - 1; ey++) {
			for (int ez = 0; ez < nz - 1; ez++) {

				double strain_aux[6] = { 0.0 };
				double stress_aux[6] = { 0.0 };

				for (int gp = 0; gp < npe; gp++) {

					double stress_gp[6], strain_gp[6];

					get_strain(gp, strain_gp, ex, ey, ez);
					get_stress(gp, strain_gp, &non_linear_flag, stress_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; v++) {
						strain_aux[v] += strain_gp[v] * wg;
						stress_aux[v] += stress_gp[v] * wg;
					}

				}

				const int e = glo_elem3D(ex, ey, ez);
				for (int v = 0; v < nvoi; v++) {
					elem_strain[e * nvoi + v] = strain_aux[v] * ivol;
					elem_stress[e * nvoi + v] = stress_aux[v] * ivol;
				}
			}
		}
	}
}

// Explicit instantiation
template class micropp<3>;
