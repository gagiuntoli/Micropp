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


#include <iostream>
#include "ell.hpp"
#include "instrument.hpp"

using namespace std;

#define cudaCheckError() { \
cudaError_t e=cudaGetLastError(); \
if(e != cudaSuccess) { \
	printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
	exit(0); \
} \
}

__global__
void ell_mvp_kernel(const ell_matrix *m_d, const double *vals_d,
		    const int *cols_d, const double *x_d, double *y_d)
{
	__shared__ double cache[2 * 128];

	const unsigned int col = threadIdx.y; // 0 .. 127
	const unsigned int row_begin = threadIdx.x + blockDim.x * blockIdx.x;
	const unsigned int stride = gridDim.x * blockDim.x;

	const unsigned int num_strides = m_d->nrow / stride + 1;

	for (int row = row_begin; row < stride * num_strides; row += stride) {

	const unsigned int ix = row * m_d->nnz + threadIdx.y;

	cache[blockDim.y * threadIdx.x + threadIdx.y] = 
		(row < m_d->nrow && col < m_d->nnz) ? vals_d[ix] * x_d[cols_d[ix]] : 0;
	__syncthreads();

	for (unsigned int s = 1; s < blockDim.y; s *= 2) {
		if (threadIdx.y % (2 * s) == 0 && threadIdx.y < m_d->nnz) {
			cache[blockDim.y * threadIdx.x + threadIdx.y] += 
			cache[blockDim.y * threadIdx.x + threadIdx.y + s];
		}
		__syncthreads();
	}
	// save result for this block on global memory
	if (threadIdx.y == 0 && row < m_d->nrow) y_d[row] = cache[blockDim.y * threadIdx.x];
        }

}

__global__
void get_diag_kernel(const ell_matrix *m_d, const double *vals_d, double *k_d)
{
	int it = threadIdx.x + blockDim.x * blockIdx.x;
	int stride_x = blockDim.x * gridDim.x;

	for (int i = it; i < m_d->nn; i += stride_x) {
		for (int d = 0; d < m_d->nfield; ++d) {
			k_d[i * m_d->nfield + d] = 1 / vals_d[i * m_d->nfield * m_d->nnz
				+ m_d->shift * m_d->nfield + d * m_d->nnz + d];
		}
	}
}

__global__
void axpby(const double alpha, const double *x, const double beta, double *y, const int n)
{
	// y = alpha x + beta y
	int it = threadIdx.x + blockDim.x * blockIdx.x;
	int stride_x = blockDim.x * gridDim.x;

	for (int i = it; i < n; i += stride_x) {
		y[i] = alpha * x[i] + beta * y[i];
	}
}

__global__
void dot_i(double *z, const double *x, const double *y, const int n)
{
	int it = threadIdx.x + blockDim.x * blockIdx.x;
	int stride_x = blockDim.x * gridDim.x;

	for (int i = it; i < n; i += stride_x) {
		z[i] = x[i] * y[i];
	}
}


__global__
void dot_prod_kernel(const double *arr1_d, const double *arr2_d, double *g_res, const int n)
{    
	// There is a different shared memory for each block
	extern __shared__ double sdata[]; 
	int it = threadIdx.x + blockDim.x * blockIdx.x;

	sdata[threadIdx.x] = (it < n) ? arr1_d[it] * arr2_d[it] : 0; // mv data to shared memory
	__syncthreads();

	for (unsigned int s = 1; s < blockDim.x; s *= 2) {
		if (threadIdx.x % (2 * s) == 0) {
			sdata[threadIdx.x] += sdata[threadIdx.x + s];
		}
		__syncthreads();
	}
	//save result for this block on global memory
	if (threadIdx.x == 0) g_res[blockIdx.x] = sdata[0];
}


int ell_solve_cgpd_cuda(const ell_matrix *m, const double *b, double *x, double *err)
{
	INST_START;

	const int grid = 512;
	const int block = 1024;

	const int grid_dot = 100000;
	const int block_dot = 512;

	ell_matrix *m_d;
	double *vals_d;
	int *cols_d;
	double *dot_d;
	double dot_h;
	double *b_d, *k_d, *r_d, *x_d, *z_d, *p_d, *Ap_d;
	double *g_res_h, *g_res_d;

	cudaMalloc((void **)&b_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&k_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&r_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&x_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&z_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&p_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&Ap_d, m->nrow * sizeof(double));

	g_res_h = (double *)malloc(grid_dot * sizeof(double));
	cudaMalloc((void **)&g_res_d, grid_dot * sizeof(double));

	cudaMalloc((void **)&m_d, sizeof(ell_matrix));
	cudaMalloc((void **)&vals_d, m->nrow * m->nnz * sizeof(double));
	cudaMalloc((void **)&cols_d, m->nrow * m->nnz * sizeof(int));
	cudaMalloc((void **)&dot_d, sizeof(double));


	cudaMemcpy(m_d, m, sizeof(ell_matrix), cudaMemcpyHostToDevice);
	cudaMemcpy(vals_d, m->vals, m->nrow * m->nnz * sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(cols_d, m->cols, m->nrow * m->nnz * sizeof(int), cudaMemcpyHostToDevice);

	/* Conjugate Gradient Algorithm (CG) with Jacobi Preconditioner */

	if (!m || !b || !x)
		return 1;

	for (int i = 0; i < m->nn; i++) {
		for (int d = 0; d < m->nfield; d++)
			m->k[i * m->nfield + d] = 1 / m->vals[i * m->nfield * m->nnz
				+ m->shift * m->nfield + d * m->nnz + d];
	}

	cudaMemset(x_d, 0, m->nrow * sizeof(double));

	cudaMemcpy(r_d, b, m->nrow * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(k_d, m->k, m->nrow * sizeof(double), cudaMemcpyHostToDevice);
	dot_i<<<grid, block>>>(z_d, k_d, r_d, m->nrow);

	cudaMemcpy(p_d, z_d, m->nrow * sizeof(double), cudaMemcpyDeviceToDevice);

	dot_prod_kernel<<<grid_dot, block_dot, block_dot* 8>>>(r_d, z_d, g_res_d, m->nrow);
	cudaMemcpy(g_res_h, g_res_d, grid_dot * sizeof(double), cudaMemcpyDeviceToHost);
	dot_h = 0.0;
	for (int i = 0; i < grid_dot; ++i) {
		dot_h += g_res_h[i];
	}
	double rz = dot_h;

	dot_prod_kernel<<<grid_dot, block_dot, block_dot* 8>>>(z_d, z_d, g_res_d, m->nrow);
	cudaMemcpy(g_res_h, g_res_d, grid_dot * sizeof(double), cudaMemcpyDeviceToHost);
	dot_h = 0.0;
	for (int i = 0; i < grid_dot; ++i) {
		dot_h += g_res_h[i];
	}

	double pnorm_0 = sqrt(dot_h);
	double pnorm = pnorm_0;
	cout << m->nrow << endl;

	int its = 0;
	while (its < m->max_its) {

		if (pnorm < m->min_err || pnorm < pnorm_0 * m->rel_err)
			break;

		dim3 grid_mvp(10000, 1, 1);
		dim3 block_mvp(2, 128, 1);
		ell_mvp_kernel<<<grid_mvp, block_mvp>>>(m_d, vals_d, cols_d, p_d, Ap_d);

		dot_prod_kernel<<<grid_dot, block_dot, block_dot* 8>>>(p_d, Ap_d, g_res_d, m->nrow);
		cudaMemcpy(g_res_h, g_res_d, grid_dot * sizeof(double), cudaMemcpyDeviceToHost);
		dot_h = 0.0;
		for (int i = 0; i < grid_dot; ++i) {
			dot_h += g_res_h[i];
		}

		double pAp = dot_h;

		const double alpha = rz / pAp;

	        axpby<<<grid, block>>>(alpha, p_d, 1.0, x_d, m->nrow);

	        axpby<<<grid, block>>>(-1.0 * alpha, Ap_d, 1.0, r_d, m->nrow);

		dot_i<<<grid, block>>>(z_d, k_d, r_d, m->nrow);

		dot_prod_kernel<<<grid_dot, block_dot, block_dot* 8>>>(z_d, z_d, g_res_d, m->nrow);
		cudaMemcpy(g_res_h, g_res_d, grid_dot * sizeof(double), cudaMemcpyDeviceToHost);
		dot_h = 0.0;
		for (int i = 0; i < grid_dot; ++i) {
			dot_h += g_res_h[i];
		}
		pnorm = sqrt(dot_h);

		dot_prod_kernel<<<grid_dot, block_dot, block_dot* 8>>>(r_d, z_d, g_res_d, m->nrow);
		cudaMemcpy(g_res_h, g_res_d, grid_dot * sizeof(double), cudaMemcpyDeviceToHost);
		dot_h = 0.0;
		for (int i = 0; i < grid_dot; ++i) {
			dot_h += g_res_h[i];
		}
		double rz_n = dot_h;

		const double beta = rz_n / rz;

	        axpby<<<grid, block>>>(1.0, z_d, beta, p_d, m->nrow);

		rz = rz_n;
		its++;
	}

	*err = rz;

	cudaMemcpy(x, x_d, m->nrow * sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(b_d);
	cudaFree(k_d);
	cudaFree(r_d);
	cudaFree(x_d);
	cudaFree(z_d);
	cudaFree(p_d);
	cudaFree(Ap_d);

	cudaFree(m_d);
	cudaFree(vals_d);
	cudaFree(cols_d);
	free(g_res_h);
	cudaFree(g_res_d);

	return its;
}
