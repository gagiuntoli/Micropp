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
	int it = threadIdx.x + blockDim.x * blockIdx.x;
	int stride_x = blockDim.x * gridDim.x;

	for (int i = it; i < m_d->nrow; i += stride_x) {
		double tmp = 0;
		const int ix = i * m_d->nnz;
		for (int j = 0; j < m_d->nnz; j++) {
			tmp += vals_d[ix + j] * x_d[cols_d[ix + j]];
		}
		y_d[i] = tmp;
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
void cpxy(double *y, const double *x, const int n)
{
	// y = x 
	int it = threadIdx.x + blockDim.x * blockIdx.x;
	int stride_x = blockDim.x * gridDim.x;

	for (int i = it; i < n; i += stride_x) {
		y[i] = x[i];
	}
}


// atomicAdd for double
__device__ double atomicAdd1_d(double* address, double val)
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


#define THREADS_PER_BLOCK 512

__global__ void dot_kernel(double *result, const double *x, const double *y, const int n)
{
	__shared__ double temp[THREADS_PER_BLOCK];

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < n) {
		temp[threadIdx.x] = x[i] * y[i];
	}

	__syncthreads();

	if (threadIdx.x == 0)
	{
		double sum = 0;
		for (int j = 0; j < THREADS_PER_BLOCK; j++) {
			sum += temp[j];
		}
		atomicAdd1_d(result, sum);
	}
}

__global__ void dot_kernel_seq(double *result, const double *x, const double *y, const int n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i == 0)
	{
		double sum = 0;
		for (int j = 0; j < n; j++) {
			sum += x[j] * y[j];
		}
		*result = sum;
	}
}


int ell_solve_cgpd_cuda_1(const ell_matrix *m_h, const double *b_h, double *x_h, double *err)
{
	INST_START;

	/* Conjugate Gradient Algorithm (CG) with Jacobi Preconditioner */

	if (!m_h || !b_h || !x_h)
		return 1;

	ell_matrix *m_d;
	double *vals_d;
	int *cols_d;
	const int n = m_h->nrow;

	double dot_h, *dot_d;
	double *b_d, *k_d, *r_d, *x_d, *z_d, *p_d, *Ap_d;

	cudaMalloc((void **)&m_d, sizeof(ell_matrix));
	cudaMalloc((void **)&vals_d, m_h->nrow * m_h->nnz * sizeof(double));
	cudaMalloc((void **)&cols_d, m_h->nrow * m_h->nnz * sizeof(int));

	cudaMalloc((void **)&dot_d, sizeof(double));
	cudaMalloc((void **)&b_d, n * sizeof(double));
	cudaMalloc((void **)&k_d, n * sizeof(double));
	cudaMalloc((void **)&r_d, n * sizeof(double));
	cudaMalloc((void **)&x_d, n * sizeof(double));
	cudaMalloc((void **)&z_d, n * sizeof(double));
	cudaMalloc((void **)&p_d, n * sizeof(double));
	cudaMalloc((void **)&Ap_d, n * sizeof(double));

	cudaMemcpy(m_d, m_h, sizeof(ell_matrix), cudaMemcpyHostToDevice);
	cudaMemcpy(vals_d, m_h->vals, m_h->nrow * m_h->nnz * sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(cols_d, m_h->cols, m_h->nrow * m_h->nnz * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(b_d, b_h, n * sizeof(double), cudaMemcpyHostToDevice);

	const int grid = 128;
	const int block = 128;

	get_diag_kernel<<<grid, block>>>(m_d, vals_d, k_d);
	//for (int i = 0; i < m->nn; i++) {
	//	for (int d = 0; d < m->nfield; d++)
	//		m->k[i * m->nfield + d] = 1 / m->vals[i * m->nfield * m->nnz
	//			+ m->shift * m->nfield + d * m->nnz + d];
	//}
	cudaMemcpy(m_h->k, k_d, n * sizeof(double), cudaMemcpyDeviceToHost);

	cudaMemset(x_d, 0, n * sizeof(double));
	//for (int i = 0; i < m->nrow; ++i)
	//	x[i] = 0.0;

	ell_mvp_kernel<<<grid, block>>>(m_d, vals_d, cols_d, x_d, r_d);
	//ell_mvp(m, x, m->r);

	axpby<<<grid, block>>>(1.0, b_d, -1.0, r_d, n);
	//for (int i = 0; i < m->nrow; ++i)
	//	m->r[i] = b[i] - m->r[i];

	dot_i<<<grid, block>>>(z_d, k_d, r_d, n);
	//for (int i = 0; i < m->nrow; ++i)
	//	m->z[i] = m->k[i] * m->r[i];
	cudaMemcpy(m_h->z, z_d, n * sizeof(double), cudaMemcpyDeviceToHost);

	cpxy<<<grid, block>>>(p_d, z_d, n);
	cudaMemcpy(m_h->p, p_d, n * sizeof(double), cudaMemcpyDeviceToHost);
	//for (int i = 0; i < m->nrow; ++i)
	//	m->p[i] = m->z[i];

	dot_kernel<<<grid, block>>>(dot_d, r_d, z_d, n);
	cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
	double rz = dot_h;
	//double rz = get_dot(m->r, m->z, m->nrow);

	dot_kernel<<<grid, block>>>(dot_d, z_d, z_d, n);
	cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
	double pnorm_0 = sqrt(dot_h);
	//double pnorm_0 = sqrt(get_dot(m->z, m->z, m->nrow));

	double pnorm = pnorm_0;
	cout << "err : " << pnorm << endl;

	int its = 0;
	while (its < m_h->max_its) {

		if (pnorm < m_h->min_err || pnorm < pnorm_0 * m_h->rel_err)
			break;

		ell_mvp_kernel<<<grid, block>>>(m_d, vals_d, cols_d, p_d, Ap_d);
		cudaMemcpy(m_h->Ap, Ap_d, n * sizeof(double), cudaMemcpyDeviceToHost);
		//ell_mvp(m, m->p, m->Ap);

		dot_kernel<<<grid, block>>>(dot_d, p_d, Ap_d, n);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
		double pAp = dot_h;
		//double pAp = get_dot(m->p, m->Ap, m->nrow);

		const double alpha = rz / pAp;

	        axpby<<<grid, block>>>(alpha, p_d, 0.0, x_d, n);
		cudaMemcpy(x_h, x_d, n * sizeof(double), cudaMemcpyDeviceToHost);
		//for (int i = 0; i < m->nrow; ++i)
		//	x[i] += alpha * m->p[i];

	        axpby<<<grid, block>>>(-1.0 * alpha, Ap_d, 0.0, r_d, n);
		cudaMemcpy(m_h->r, r_d, n * sizeof(double), cudaMemcpyDeviceToHost);
		//for (int i = 0; i < m->nrow; ++i)
		//	m->r[i] -= alpha * m->Ap[i];

		dot_i<<<grid, block>>>(z_d, k_d, r_d, n);
		cudaMemcpy(m_h->z, z_d, n * sizeof(double), cudaMemcpyDeviceToHost);
		//for (int i = 0; i < m->nrow; ++i)
		//	m->z[i] = m->k[i] * m->r[i];

		dot_kernel<<<grid, block>>>(dot_d, z_d, z_d, n);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
		pnorm = sqrt(dot_h);
		//pnorm = sqrt(get_dot(m->z, m->z, m->nrow));

		dot_kernel<<<grid, block>>>(dot_d, r_d, z_d, n);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
		double rz_n = dot_h;
		//double rz_n = get_dot(m->r, m->z, m->nrow);

		const double beta = rz_n / rz;
	        axpby<<<grid, block>>>(1.0, z_d, beta, p_d, n);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
		//for (int i = 0; i < m->nrow; ++i)
		//	m->p[i] = m->z[i] + beta * m->p[i];

		cout << "err : " << rz << endl;
		rz = rz_n;
		its++;
	}

	*err = rz;

	cudaMemcpy(x_h, x_d, n * sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(m_d);
	cudaFree(vals_d);
	cudaFree(cols_d);

	cudaFree(dot_d);
	cudaFree(b_d);
	cudaFree(k_d);
	cudaFree(r_d);
	cudaFree(x_d);
	cudaFree(z_d);
	cudaFree(p_d);
	cudaFree(Ap_d);

	return its;
}

int ell_solve_cgpd_cuda(const ell_matrix *m, const double *b, double *x, double *err)
{
	INST_START;

	const int grid = 256;
	const int block = 512;

	ell_matrix *m_d;
	double *vals_d;
	int *cols_d;
	double *dot_d;
	double dot_h;
	double *b_d, *k_d, *r_d, *x_d, *z_d, *p_d, *Ap_d;

	cudaMalloc((void **)&b_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&k_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&r_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&x_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&z_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&p_d, m->nrow * sizeof(double));
	cudaMalloc((void **)&Ap_d, m->nrow * sizeof(double));

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
	for (int i = 0; i < m->nrow; ++i)
		x[i] = 0.0;

	ell_mvp(m, x, m->r);

	for (int i = 0; i < m->nrow; ++i)
		m->r[i] = b[i] - m->r[i];

	for (int i = 0; i < m->nrow; ++i)
		m->z[i] = m->k[i] * m->r[i];

	for (int i = 0; i < m->nrow; ++i)
		m->p[i] = m->z[i];

	cudaMemcpy(k_d, m->k, m->nrow * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(r_d, m->r, m->nrow * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(z_d, m->z, m->nrow * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(p_d, m->p, m->nrow * sizeof(double), cudaMemcpyHostToDevice);

	double rz = get_dot(m->r, m->z, m->nrow);

	double pnorm_0 = sqrt(get_dot(m->z, m->z, m->nrow));
	double pnorm = pnorm_0;

	int its = 0;
	while (its < m->max_its) {

		if (pnorm < m->min_err || pnorm < pnorm_0 * m->rel_err)
			break;

		ell_mvp_kernel<<<grid, block>>>(m_d, vals_d, cols_d, p_d, Ap_d);

//dot_kernel<<<m->nrow / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK>>>(dot_d, p_d, Ap_d, m->nrow);
		dot_kernel_seq<<<1, 1>>>(dot_d, p_d, Ap_d, m->nrow);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
		double pAp = dot_h;

		const double alpha = rz / pAp;

	        axpby<<<grid, block>>>(alpha, p_d, 1.0, x_d, m->nrow);

	        axpby<<<grid, block>>>(-1.0 * alpha, Ap_d, 1.0, r_d, m->nrow);

		dot_i<<<grid, block>>>(z_d, k_d, r_d, m->nrow);

		dot_kernel_seq<<<1, 1>>>(dot_d, z_d, z_d, m->nrow);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
		pnorm = sqrt(dot_h);

		dot_kernel_seq<<<1, 1>>>(dot_d, r_d, z_d, m->nrow);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
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

	return its;
}
