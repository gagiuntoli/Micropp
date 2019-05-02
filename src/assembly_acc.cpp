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
void micropp<3>::get_stress_acc(int gp, const double eps[nvoi],
				const double *vars_old,
				double stress_gp[nvoi],
				int ex, int ey, int ez) const
{
	const int e = glo_elem(ex, ey, ez);
	const material_t *material = get_material(e);
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
	const material_t *material = get_material(e);

	double ctan[nvoi][nvoi];
	constexpr int npedim = npe * dim;
	constexpr int npedim2 = npedim * npedim;

	double TAe[npedim2] = { 0.0 };

	for (int gp = 0; gp < npe; ++gp) {

		double eps[6];
		get_strain(u, gp, eps, ex, ey, ez);

		const double *vars = (vars_old) ? &vars_old[intvar_ix(e, gp, 0)] : nullptr;
		material->get_ctan(eps, (double *)ctan, vars);

		double cxb[nvoi][npedim];

		for (int i = 0; i < nvoi; ++i) {
			for (int j = 0; j < npedim; ++j) {
				double tmp = 0.0;
				for (int k = 0; k < nvoi; ++k)
					tmp += ctan[i][k] * calc_bmat_cache[gp][k][j];
				cxb[i][j] = tmp * wg;
			}
		}

		for (int m = 0; m < nvoi; ++m) {
			for (int i = 0; i < npedim; ++i) {
				const int inpedim = i * npedim;
				const double bmatmi = calc_bmat_cache[gp][m][i];
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
	double stress_gp[nvoi], strain_gp[nvoi];

	memset(be, 0, npedim * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		get_strain(u, gp, strain_gp, ex, ey, ez);
		get_stress_acc(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);

		for (int i = 0; i < npedim; ++i)
			for (int j = 0; j < nvoi; ++j)
				be[i] += calc_bmat_cache[gp][j][i] * stress_gp[j] * wg;
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

	double* strain_gp = new double [nex*ney*nez*npe*nvoi];

#pragma acc parallel loop copy(strain_gp[:nex*ney*nez*npe*nvoi]) copyin(u[:nndim])
	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {
				for (int gp = 0; gp < npe; ++gp) {
					get_strain(u, gp, &strain_gp[ex*ney*nez*npe*nvoi+ey*nez*npe*nvoi+ez*npe*nvoi+gp*nvoi], ex, ey, ez);
				}
			}
		}
	}
	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				int n[npe];
				get_elem_nodes(n, ex, ey, ez);

				for (int j = 0; j < npe; ++j)
					for (int d = 0; d < dim; ++d)
						index[j * dim + d] = n[j] * dim + d;

				constexpr int npedim = npe * dim;
				double stress_gp[nvoi];

				memset(be, 0, npedim * sizeof(double));

				for (int gp = 0; gp < npe; ++gp) {

					get_stress_acc(gp, &strain_gp[ex*ney*nez*npe*nvoi+ey*nez*npe*nvoi+ez*npe*nvoi+gp*nvoi], vars_old, stress_gp, ex, ey, ez);

					for (int i = 0; i < npedim; ++i)
						for (int j = 0; j < nvoi; ++j)
							be[i] += calc_bmat_cache[gp][j][i] * stress_gp[j] * wg;
				}

				for (int i = 0; i < npe * dim; ++i)
					b[index[i]] += be[i];
			}
		}
	}
	delete[] strain_gp;

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

	double *bmat = new double[npe*nvoi*npedim];
	for (int gp = 0; gp < npe; ++gp) {
		for (int i = 0; i < nvoi; i++){
		  for (int j = 0; j < npedim; j++){
				bmat[gp*nvoi*npedim+i*npedim+j] = calc_bmat_cache[gp][i][j];
	    }
    }
	}
	const int nnz_nfield = nfield * nnz;
	const int npe_nfield = npe * nfield;
	const int npe_nfield2 = npe * nfield * nfield;

	double* eps = new double[nex*ney*nez*npe*6];
#pragma acc parallel loop copy(eps[:nex*ney*nez*npe*6]) copyin(u[:nndim])
	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {
				for (int gp = 0; gp < npe; ++gp) {
					get_strain(u, gp, &eps[ex*ney*nez*npe*6+ey*nez*npe*6+ez*npe*6+gp*6], ex, ey, ez);
				}
			}
		}
	}
	double *ctan = new double[nex*ney*nez*npe*nvoi*nvoi];
	int *ix_glo = new int[nex*ney*nez*8];
	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {
				const int e = glo_elem(ex, ey, ez);
				const material_t *material = get_material(e);
				for (int gp = 0; gp < npe; ++gp) {
					const double *vars = (vars_old) ? &vars_old[intvar_ix(e, gp, 0)] : nullptr;
					material->get_ctan(&eps[ex*ney*nez*npe*6+ey*nez*npe*6+ez*npe*6+gp*6], &ctan[ex*ney*nez*npe*nvoi*nvoi+ey*nez*npe*nvoi*nvoi+ez*npe*nvoi*nvoi+gp*nvoi*nvoi], vars);
				}
				const int n0 = ez * nxny + ey * nx + ex;
				const int n1 = n0 + 1;
				const int n2 = n0 + nx + 1;
				const int n3 = n0 + nx;

				ix_glo[ex*ney*nez*8+ey*nez*8+ez*8+0] =n0;
				ix_glo[ex*ney*nez*8+ey*nez*8+ez*8+1] =n1;
				ix_glo[ex*ney*nez*8+ey*nez*8+ez*8+2] =n2;
				ix_glo[ex*ney*nez*8+ey*nez*8+ez*8+3] =n3;
				ix_glo[ex*ney*nez*8+ey*nez*8+ez*8+4] =n0 + nxny;
				ix_glo[ex*ney*nez*8+ey*nez*8+ez*8+5] =n1 + nxny;
				ix_glo[ex*ney*nez*8+ey*nez*8+ez*8+6] =n2 + nxny;
				ix_glo[ex*ney*nez*8+ey*nez*8+ez*8+7] =n3 + nxny;
			}
		}
	}
	delete[] eps;
#pragma acc parallel loop gang vector copyin(ctan[:nex*ney*nez*npe*nvoi*nvoi],bmat[:npe*nvoi*npedim],A[:1],A->nrow,A->nnz,cols_row[:8][:8],ix_glo[:nex*ney*nez*8])copy(A->vals[:A->nrow * A->nnz])
	for (int ex = 0; ex < nex*ney*nez; ++ex) {

				double Ae[npedim2];
				for(int i=0;i<npedim2;i++)Ae[i]=0;

				for (int gp = 0; gp < npe; ++gp) {
					double cxb[nvoi][npedim];
					for (int i = 0; i < nvoi; ++i) {
						for (int j = 0; j < npedim; ++j) {
							double tmp = 0.0;
							for (int k = 0; k < nvoi; ++k){
								tmp += ctan[ex*npe*nvoi*nvoi+gp*nvoi*nvoi+i*nvoi+k] * bmat[gp*nvoi*npedim+k*npedim+j];
							}
							cxb[i][j] = tmp * wg;
						}
					}
					for (int m = 0; m < nvoi; ++m) {
						for (int i = 0; i < npedim; ++i) {
							const int inpedim = i * npedim;
							const double bmatmi = bmat[gp*nvoi*npedim+m*npedim+i];
							for (int j = 0; j < npedim; ++j){
								Ae[inpedim + j] += bmatmi * cxb[m][j];
							}
						}
					}
				}

				for (int fi = 0; fi < nfield; ++fi){
					for (int fj = 0; fj < nfield; ++fj){
						for (int i = 0; i < npe; ++i){
							for (int j = 0; j < npe; ++j){
#pragma acc atomic update
								A->vals[ix_glo[ex*8+i] * nnz_nfield + cols_row[i][j] * nfield + fi * nnz + fj] += Ae[i * npe_nfield2 + fi * npe_nfield + j * nfield + fj];
							}
						}
					}
				}
	}
	delete []ix_glo;
	delete []bmat;
	delete []ctan;
	ell_set_bc_3D_acc(A);
}


