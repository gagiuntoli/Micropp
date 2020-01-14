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


#include <cstdio>
#include "micropp.hpp"
#include "ell.hpp"
#include "instrument.hpp"
#include "common.hpp"

#include "cuda_profiler_api.h"

#define cudaCheckError() { \
cudaError_t e=cudaGetLastError(); \
if(e != cudaSuccess) { \
	printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
	exit(0); \
} \
}

#define DIM 3
#define NPE 8
#define NVOI 6
#define NVOI2 (NVOI * NVOI)
#define NPEDIM (NPE * DIM)
#define NPEDIM2 (NPEDIM * NPEDIM)

using namespace std;

// Structute to pass parameters to GPU
struct Params {

	int dim, npe, nvoi;
	int nex, ney, nez;
	int nx, ny, nz;
	double bmat_cache[NPE][NVOI][NPE * DIM];
	double wg;

};

// atomicAdd for double
__device__ double atomicAdd_d(double* address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do{
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, 
				__double_as_longlong(val + __longlong_as_double(assumed)));
	} while(assumed != old);
	return __longlong_as_double(old);
}


__device__
void ell_add_3D_gpu(ell_matrix *m, double *vals_d, int ex, int ey, int ez, 
		    const double Ae[8 * 3 * 8 * 3])
{
	// assembly Ae in 3D structured grid representation
	// nFields : number of scalar components on each node

	const int nx = m->n[0];
	const int ny = m->n[1];
	const int nfield = m->nfield;
	const int npe = 8;
	const int nnz = m->nnz;
	const int cols_row[8][8] = {
		{ 13, 14, 17, 16, 22, 23, 26, 25 },
		{ 12, 13, 16, 15, 21, 22, 25, 24 },
		{ 9,  10, 13, 12, 18, 19, 22, 21 },
		{ 10, 11, 14, 13, 19, 20, 23, 22 },
		{ 4,  5,  8,  7,  13, 14, 17, 16 },
		{ 3,  4,  7,  6,  12, 13, 16, 15 },
		{ 0,  1,  4,  3,  9,  10, 13, 12 },
		{ 1,  2,  5,  4,  10, 11, 14, 13 } };

	const int nxny = nx * ny;
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

	for (int fi = 0; fi < nfield; ++fi)
		for (int fj = 0; fj < nfield; ++fj)
			for (int i = 0; i < npe; ++i)
				for (int j = 0; j < npe; ++j){
					atomicAdd_d(&vals_d[ix_glo[i] * nnz_nfield + cols_row[i][j] * nfield + fi * nnz + fj],
						    Ae[i * npe_nfield2 + fi * npe_nfield + j * nfield + fj]);
				}

}

__device__
void get_ctan_d(const double *eps, double *ctan, const double *history_params)
{
	const double E = 1.0e7;
	const double nu = 0.25;

	const double lambda = nu * E / ((1. + nu) * (1. - 2. * nu));
	const double mu = E / (2. * (1. + nu));

	memset(ctan, 0, 6 * 6 * sizeof(double));

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			ctan[i * 6 + j] += lambda;

	for (int i = 0; i < 3; ++i)
		ctan[i * 6 + i] += 2 * mu;

	for (int i = 3; i < 6; ++i)
		ctan[i * 6 + i] = mu;
}

__device__
void get_stress_d(const double eps[NVOI], double stress[NVOI], const double *history_params)
{
	const double E = 1.0e7;
	const double nu = 0.25;

	const double lambda = nu * E / ((1. + nu) * (1. - 2. * nu));
	const double mu = E / (2. * (1. + nu));

	// stress[i][j] = lambda eps[k][k] * delta[i][j] + mu eps[i][j]
	for (int i = 0; i < 3; ++i)
		stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
			    + 2 * mu * eps[i];

	for (int i = 3; i < 6; ++i)
		stress[i] = mu * eps[i];
}


__global__
void assembly_kernel(ell_matrix *A_d, double *vals_d, const double *u, Params *params_d)
{
	const int nex = params_d->nex;
	const int ney = params_d->ney;
	const int nez = params_d->nez;
	const double wg = params_d->wg;

	int ex_t = threadIdx.x + blockDim.x * blockIdx.x;
	int ey_t = threadIdx.y + blockDim.y * blockIdx.y;
	int ez_t = threadIdx.z + blockDim.z * blockIdx.z;
	int stride_x = blockDim.x * gridDim.x;
	int stride_y = blockDim.y * gridDim.y;
	int stride_z = blockDim.z * gridDim.z;

	for (int ex = ex_t; ex < nex; ex += stride_x) {
		for (int ey = ey_t; ey < ney; ey += stride_y) {
			for (int ez = ez_t; ez < nez; ez += stride_z) {
	double TAe[NPEDIM2] = { 0.0 };
	for (int gp = 0; gp < NPE; ++gp) {

		double eps[NVOI];
		double ctan[NVOI2];
		get_strain(u, gp, eps, params_d->bmat_cache, params_d->nx, params_d->ny, ex, ey, ez);
		get_ctan_d(eps, ctan, nullptr);
		double cxb[NVOI][NPEDIM];

		for (int i = 0; i < NVOI; ++i) {
			for (int j = 0; j < NPEDIM; ++j) {
				double tmp = 0.0;
				for (int k = 0; k < NVOI; ++k)
					tmp += ctan[i * NVOI + k] 
						* params_d->bmat_cache[gp][k][j];
				cxb[i][j] = tmp * wg;
			}
		}

		for (int m = 0; m < NVOI; ++m) {
			for (int i = 0; i < NPEDIM; ++i) {
				const int inpedim = i * NPEDIM;
				const double bmatmi = params_d->bmat_cache[gp][m][i];
				for (int j = 0; j < NPEDIM; ++j)
					TAe[inpedim + j] += bmatmi * cxb[m][j];
			}
		}
	}
	ell_add_3D_gpu(A_d, vals_d, ex, ey, ez, TAe);
			}
		}
	}
}


template <>
void micropp<3>::assembly_mat(ell_matrix *A, const double *u, const double *vars_old)
{
	INST_START;

	cout << "assembly_mat_cuda" << endl;
	//cudaProfilerStart();
	ell_set_zero_mat(A);

	Params params_h;

	params_h.nx = nx;
	params_h.ny = ny;
	params_h.nz = nz;
	params_h.nex = nex;
	params_h.ney = ney;
	params_h.nez = nez;
	params_h.wg = wg;
	memcpy(&params_h.bmat_cache, bmat_cache, NPE * NVOI * NPE * DIM * sizeof(double));

	Params *params_d;
	double *u_d;
	double *vals_d;
	ell_matrix *A_d;

	cudaMalloc((void **)&params_d, sizeof(Params));
	cudaMalloc((void**)&A_d, sizeof(ell_matrix));
	cudaMalloc((void**)&vals_d, A->nnz * A->nrow * sizeof(double));
	cudaMalloc((void**)&u_d, nndim * sizeof(double));
	cudaMemcpy(params_d, &params_h, sizeof(Params), cudaMemcpyHostToDevice);
	cudaMemcpy(u_d, u, nndim * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(A_d, A, sizeof(ell_matrix), cudaMemcpyHostToDevice);
	cudaMemset(vals_d, 0, A->nnz * A->nrow * sizeof(double));

	dim3 grid(15, 15, 15);
	dim3 block(4, 4, 4);
	assembly_kernel<<<grid, block>>>(A_d, vals_d, u_d, params_d);
        cudaCheckError();

	cudaMemcpy(A->vals, vals_d, A->nrow * A->nnz * sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(A_d);
	cudaFree(vals_d);
	cudaFree(u_d);
	cudaFree(params_d);

	ell_set_bc_3D(A);
	//cudaProfilerStop();

	return;
}

__device__
void get_elem_rhs(const double *u, const double *vars_old, double be[NPEDIM],
		  Params *params_d, int ex, int ey, int ez)
{
	double stress_gp[NVOI], strain_gp[NVOI];

	memset(be, 0, NPEDIM * sizeof(double));

	for (int gp = 0; gp < NPE; ++gp) {

		get_strain(u, gp, strain_gp, params_d->bmat_cache, params_d->nx, params_d->ny, ex, ey, ez);
		get_stress_d(strain_gp, stress_gp, vars_old);

		for (int i = 0; i < NPEDIM; ++i)
			for (int j = 0; j < NVOI; ++j)
				be[i] += params_d->bmat_cache[gp][j][i] 
					* stress_gp[j] * params_d->wg;
	}
}

__global__
void assembly_rhs_kernel(double *b_d, const double *u, Params *params_d)
{
	const int nex = params_d->nex;
	const int ney = params_d->ney;
	const int nez = params_d->nez;

	int ex_t = threadIdx.x + blockDim.x * blockIdx.x;
	int ey_t = threadIdx.y + blockDim.y * blockIdx.y;
	int ez_t = threadIdx.z + blockDim.z * blockIdx.z;
	int stride_x = blockDim.x * gridDim.x;
	int stride_y = blockDim.y * gridDim.y;
	int stride_z = blockDim.z * gridDim.z;

	for (int ex = ex_t; ex < nex; ex += stride_x) {
	for (int ey = ey_t; ey < ney; ey += stride_y) {
	for (int ez = ez_t; ez < nez; ez += stride_z) {

		int n[NPE];
		get_elem_nodes(n, params_d->nx, params_d->ny, ex, ey, ez);

		double be[DIM * NPE];
		int index[DIM * NPE];
		for (int j = 0; j < NPE; ++j)
			for (int d = 0; d < DIM; ++d)
				index[j * DIM + d] = n[j] * DIM + d;

		for (int gp = 0; gp < NPE; ++gp) {
			get_elem_rhs(u, nullptr, be, params_d, ex, ey, ez);
		}

		for (int i = 0; i < NPE * DIM; ++i)
			atomicAdd_d(&b_d[index[i]], be[i]);

	}
	}
	}
}

template<>
double micropp<3>::assembly_rhs(const double *u, const double *vars_old, double *b)
{
	INST_START;

	cout << "assembly_rhs_cuda" << endl;
	Params params_h;

	params_h.nx = nx;
	params_h.ny = ny;
	params_h.nz = nz;
	params_h.nex = nex;
	params_h.ney = ney;
	params_h.nez = nez;
	params_h.wg = wg;
	memcpy(&params_h.bmat_cache, bmat_cache, NPE * NVOI * NPE * DIM * sizeof(double));

	Params *params_d;
	double *u_d;
	double *b_d;

	cudaMalloc((void **)&params_d, sizeof(Params));
	cudaMalloc((void**)&u_d, nndim * sizeof(double));
	cudaMalloc((void**)&b_d, nndim * sizeof(double));

	cudaMemcpy(params_d, &params_h, sizeof(Params), cudaMemcpyHostToDevice);
	cudaMemcpy(u_d, u, nndim * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemset(b_d, 0, nndim * sizeof(double));

	dim3 grid(15, 15, 15);
	dim3 block(4, 4, 4);
	assembly_rhs_kernel<<<grid, block>>>(b_d, u, params_d);

	cudaMemcpy(b, b_d, nndim * sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(b_d);
	cudaFree(u_d);
	cudaFree(params_d);

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
