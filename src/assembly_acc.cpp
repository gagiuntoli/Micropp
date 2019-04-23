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
material_acc *micropp<3>::get_material_acc(const int e) const
{
	return material_acc_list[elem_type[e]];
}

template <>
void micropp<3>::get_stress_acc(int gp, const double eps[nvoi],
				const double *vars_old,
				double stress_gp[nvoi],
				int ex, int ey, int ez) const
{
	const int e = glo_elem(ex, ey, ez);
	const material_acc *material = get_material_acc(e);
	const double *vars = (vars_old) ? &vars_old[intvar_ix(e, gp, 0)] : nullptr;

	material->get_stress(eps, stress_gp, vars);
}


template <>
void micropp<3>::get_elem_mat_acc(const double *u,
				  const double *vars_old,
				  double Ae[npe * dim * npe * dim],
				  int ex, int ey, int ez) const
{
	INST_START;
	const int e = glo_elem(ex, ey, ez);
	const material_acc *material = get_material_acc(e);

	double ctan[nvoi][nvoi];
	constexpr int npedim = npe * dim;
	constexpr int npedim2 = npedim * npedim;

	double TAe[npedim2] = { 0.0 };

	for (int gp = 0; gp < npe; ++gp) {

		double eps[6];
		get_strain(u, gp, eps, ex, ey, ez);

		const double *vars = (vars_old) ? &vars_old[intvar_ix(e, gp, 0)] : nullptr;
		material->get_ctan(eps, (double *)ctan, vars);

		double bmat[nvoi][npedim], cxb[nvoi][npedim];
		calc_bmat(gp, bmat);

		for (int i = 0; i < nvoi; ++i) {
			for (int j = 0; j < npedim; ++j) {
				double tmp = 0.0;
				for (int k = 0; k < nvoi; ++k)
					tmp += ctan[i][k] * bmat[k][j];
				cxb[i][j] = tmp * wg;
			}
		}

		for (int m = 0; m < nvoi; ++m) {
			for (int i = 0; i < npedim; ++i) {
				const int inpedim = i * npedim;
				const double bmatmi = bmat[m][i];
				for (int j = 0; j < npedim; ++j)
					TAe[inpedim + j] += bmatmi * cxb[m][j];
			}
		}
	}
	memcpy(Ae, TAe, npedim2 * sizeof(double));
}


template <>
void micropp<3>::get_elem_rhs_acc(const double *u, const double *vars_old,
				  double be[npe * dim],
				  int ex, int ey, int ez) const
{
	INST_START;

	constexpr int npedim = npe * dim;
	double bmat[nvoi][npedim], stress_gp[nvoi], strain_gp[nvoi];

	memset(be, 0, npedim * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		calc_bmat(gp, bmat);

		get_strain(u, gp, strain_gp, ex, ey, ez);
		get_stress_acc(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);

		for (int i = 0; i < npedim; ++i)
			for (int j = 0; j < nvoi; ++j)
				be[i] += bmat[j][i] * stress_gp[j] * wg;
	}
}


template<>
double micropp<3>::assembly_rhs_acc(const double *u, const double *vars_old,
				    double *b)
{
	INST_START;

	memset(b, 0., nndim * sizeof(double));

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

				get_elem_rhs_acc(u, vars_old, be, ex, ey, ez);

				for (int i = 0; i < npe * dim; ++i)
					b[index[i]] += be[i];
			}
		}
	}

	// boundary conditions
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, 0); // z = 0
			memset(&b[n * dim], 0., dim * sizeof(double));
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, nz - 1); // z = lx
			memset(&b[n * dim], 0., dim * sizeof(double));
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, 0, k); // y = 0
			memset(&b[n * dim], 0., dim * sizeof(double));
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, ny - 1, k); // y = ly
			memset(&b[n * dim], 0., dim * sizeof(double));
		}
	}

	for (int j = 1; j < ny - 1; ++j) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(0, j, k); // x = 0
			memset(&b[n * dim], 0., dim * sizeof(double));
		}
	}

	for (int j = 1; j < ny - 1; j++) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(nx - 1, j, k); // x = lx
			memset(&b[n * dim], 0., dim * sizeof(double));
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
void micropp<3>::assembly_mat_acc(ell_matrix *A, const double *u,
				  const double *vars_old)
{
	INST_START;

	ell_set_zero_mat(A);
	const int nx = A->n[0];
	const int ny = A->n[1];
	const int nxny = nx * ny;
	const int nfield = A->nfield;
	const int npe = 8;
	const int nnz = A->nnz;
	const int cols_row[8][8] = {
		{ 13, 14, 17, 16, 22, 23, 26, 25 },
		{ 12, 13, 16, 15, 21, 22, 25, 24 },
		{ 9,  10, 13, 12, 18, 19, 22, 21 },
		{ 10, 11, 14, 13, 19, 20, 23, 22 },
		{ 4,  5,  8,  7,  13, 14, 17, 16 },
		{ 3,  4,  7,  6,  12, 13, 16, 15 },
		{ 0,  1,  4,  3,  9,  10, 13, 12 },
		{ 1,  2,  5,  4,  10, 11, 14, 13 } };

	constexpr int npedim = npe * dim;
	constexpr int npedim2 = npedim * npedim;

	double *bmat = new double[npe * nvoi * npedim];
	for (int gp = 0; gp < npe; ++gp) {
		double bmat_temp[nvoi][npedim];
		calc_bmat(gp, bmat_temp);
		for (int i = 0; i < nvoi; i++){
			for (int j = 0; j < npedim; j++){
				bmat[gp*nvoi*npedim+i*npedim+j] = bmat_temp[i][j];
			}
		}
	}

	double *ctan = new double[nex * ney * nez * npe * nvoi * nvoi];
	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {
				const int e = glo_elem(ex, ey, ez);
				const material_acc *material = get_material_acc(e);
				for (int gp = 0; gp < npe; ++gp) {
					double eps[6];
					get_strain(u, gp, eps, ex, ey, ez);
					const double *vars = (vars_old) ? &vars_old[intvar_ix(e, gp, 0)] : nullptr;
					material->get_ctan(eps, &ctan[ex*ney*nez*npe*nvoi*nvoi+ey*nez*npe*nvoi*nvoi+ez*npe*nvoi*nvoi+gp*nvoi*nvoi], vars);
				}
			}
		}
	}

#pragma acc enter data copyin(ctan[:nex*ney*nez*npe*nvoi*nvoi], bmat[:npe*nvoi*npedim])

	double *Aes = new double[nex * ney * nez * npedim2];
	double cxb[nvoi][npedim];
#pragma acc enter data copyin(cxb[:6][:24])
#pragma acc enter data create(Aes[:nex * ney * nez * npedim2])

#pragma acc parallel loop present (Aes[:nex * ney * nez * npedim2])
	for(int i = 0; i < nex * ney * nez * npedim2; ++i)
		Aes[i] = 0.0;

	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {

				const int ctan_ix = ex*ney*nez*npe*nvoi*nvoi + \
						    ey*nez*npe*nvoi*nvoi + \
						    ez*npe*nvoi*nvoi;

				const int Aes_ix = (ez * nex * ney + ey * nex + ex) * npedim2;

				for (int gp = 0; gp < npe; ++gp) {

					const int gp1_ix = gp * nvoi * nvoi;
					const int gp2_ix = gp * nvoi * npedim;

					const double *ctan_ptr = &ctan[ctan_ix + gp1_ix];
					const double *bmat_ptr = &bmat[gp2_ix];

#pragma acc parallel loop gang vector collapse(2) present(cxb[:6][:24], ctan_ptr[:36], bmat_ptr[:6*24])
					for (int i = 0; i < nvoi; ++i) {
						for (int j = 0; j < npedim; ++j) {
							double tmp = 0.0;
#pragma acc loop reduction(+:tmp)
							for (int k = 0; k < nvoi; ++k){
								const double val1 = ctan_ptr[i * nvoi + k];
								const double val2 = bmat_ptr[k * npedim + j];
								tmp += val1 * val2;
							}
							cxb[i][j] = tmp * wg;
						}
					}

#pragma acc parallel loop gang vector collapse(2) present(cxb[:6][:24], bmat_ptr[:6*24], Aes[:nex * ney * nez * npedim2])
					for (int i = 0; i < npedim; ++i) {
						for (int j = 0; j < npedim; ++j){
							double tmp = 0.0;
#pragma acc loop reduction(+:tmp)
							for (int m = 0; m < nvoi; ++m) {
								const double val1 = bmat_ptr[m * npedim + i];
								const double val2 = cxb[m][j];
								tmp += val1 * val2;
							}
							Aes[Aes_ix + i * npedim + j] += tmp;
						}
					}
				}
			}
		}
	}

#pragma acc update self (Aes[:nex * ney * nez * npedim2])

	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {

				const int Aes_ix = (ez * nex * ney + ey * nex + ex) * npedim2;

				const int n0 = ez * nxny + ey * nx + ex;
				const int n1 = n0 + 1;
				const int n2 = n0 + nx + 1;
				const int n3 = n0 + nx;

				const int ix_glo[8] = {	n0, n1, n2, n3,
					n0 + nxny,
					n1 + nxny,
					n2 + nxny,
					n3 + nxny };

				const int nnz_nfield = nfield * nnz;
				const int npe_nfield = npe * nfield;
				const int npe_nfield2 = npe * nfield * nfield;

				for (int fi = 0; fi < nfield; ++fi){
					for (int fj = 0; fj < nfield; ++fj){
						for (int i = 0; i < npe; ++i){
							for (int j = 0; j < npe; ++j){
								A->vals[ix_glo[i] * nnz_nfield +\
								       	cols_row[i][j] * nfield + fi * nnz + fj] += \
									Aes[Aes_ix + i * npe_nfield2 + fi * npe_nfield + \
								       	j * nfield + fj];
							}
						}
					}
				}
			}
		}
	}
#pragma acc exit data delete(ctan[:nex*ney*nez*npe*nvoi*nvoi], bmat[:npe*nvoi*npedim])
#pragma acc exit data delete(cxb[:6][:24])
#pragma acc exit data delete(Aes[:nex * ney * nez * npedim2])
	delete []bmat;
	delete []ctan;
	delete []Aes;
	ell_set_bc_3D_acc(A);
}
