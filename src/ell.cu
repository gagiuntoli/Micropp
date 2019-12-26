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


#include "ell.hpp"
#include "instrument.hpp"

__global__
void ell_mvp_kernel(const ell_matrix *m_d, const double *vals_d,
		    const int *cols_d, const double *x_d, double *y_d)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x;

	if (i < m_d->nrow) {
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
	int i = threadIdx.x + blockDim.x * blockIdx.x;

	if (i < m_d->nrow) {
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
	int i = threadIdx.x + blockDim.x * blockIdx.x;

	if (i < n) {
		y[i] = alpha * x[i] + beta * y[i];
	}
}

__global__
void axpbyz(const double alpha, const double *x, 
	    const double beta, const double *y, double *z, const int n)
{
	// z = alpha x + beta y
	int i = threadIdx.x + blockDim.x * blockIdx.x;

	if (i < n) {
		z[i] = alpha * x[i] + beta * y[i];
	}
}

__global__
void dot_i(const double *x, const double *y, double *z, const int n)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x;

	if (i < n) {
		z[i] = x[i] * y[i];
	}
}

__global__
void cpxy(double *y, const double *x, const int n)
{
	// y = x 
	int i = threadIdx.x + blockDim.x * blockIdx.x;

	if (i < n) {
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


__global__ void dot_kernel(double *result, const double *x, const double *y, const int n)
{
	__shared__ int temp[512];

	*result = 0;

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < n) {
		temp[threadIdx.x] = x[i] * y[i];
	}

	__syncthreads();

	if (threadIdx.x == 0)
	{
		double sum = 0;
		for (int i = 0; i < n; i++) {
			sum += temp[i];
		}
		atomicAdd1_d(result, sum);
	}
}


int ell_solve_cgpd_cuda(const ell_matrix *m_h, const double *b_h, double *x_h, double *err)
{
	INST_START;

	ell_matrix *m_d = nullptr;
	double *vals_d = nullptr;
	int *cols_d = nullptr;
	const int n = m_h->nrow;

	double dot_h, *dot_d;
	double *b_d = nullptr;
	double *k_d, *r_d, *x_d, *z_d, *p_d, *Ap_d;

	cudaMalloc((void **)&dot_d, sizeof(double));
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

	const int grid = 64;
	const int block = 32;

	/* Conjugate Gradient Algorithm (CG) with Jacobi Preconditioner */

	if (!m_h || !b_h || !x_h)
		return 1;

	get_diag_kernel<<<grid, block>>>(m_d, vals_d, k_d);
	//for (int i = 0; i < m->nn; i++) {
	//	for (int d = 0; d < m->nfield; d++)
	//		m->k[i * m->nfield + d] = 1 / m->vals[i * m->nfield * m->nnz
	//			+ m->shift * m->nfield + d * m->nnz + d];
	//}

	cudaMemset(x_d, 0, n * sizeof(double));
	//for (int i = 0; i < m->nrow; ++i)
	//	x[i] = 0.0;

	ell_mvp_kernel<<<grid, block>>>(m_d, vals_d, cols_d, x_d, r_d);
	//ell_mvp(m, x, m->r);

	axpby<<<grid, block>>>(1.0, b_d, -1.0, r_d, n);
	//for (int i = 0; i < m->nrow; ++i)
	//	m->r[i] = b[i] - m->r[i];

	dot_i<<<grid, block>>>(k_d, r_d, z_d, n);
	//for (int i = 0; i < m->nrow; ++i)
	//	m->z[i] = m->k[i] * m->r[i];

	cpxy<<<grid, block>>>(p_d, z_d, n);
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

	int its = 0;
	while (its < m_h->max_its) {

		if (pnorm < m_h->min_err || pnorm < pnorm_0 * m_h->rel_err)
			break;

		ell_mvp_kernel<<<grid, block>>>(m_d, vals_d, cols_d, p_d, Ap_d);
		//ell_mvp(m, m->p, m->Ap);

		dot_kernel<<<grid, block>>>(dot_d, p_d, Ap_d, n);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
		double pAp = dot_h;
		//double pAp = get_dot(m->p, m->Ap, m->nrow);

		const double alpha = rz / pAp;

	        axpby<<<grid, block>>>(alpha, b_d, 0.0, x_d, n);
		//for (int i = 0; i < m->nrow; ++i)
		//	x[i] += alpha * m->p[i];

	        axpby<<<grid, block>>>(-1.0 * alpha, Ap_d, 0.0, r_d, n);
		//for (int i = 0; i < m->nrow; ++i)
		//	m->r[i] -= alpha * m->Ap[i];

		dot_i<<<grid, block>>>(k_d, r_d, z_d, n);
		//for (int i = 0; i < m->nrow; ++i)
		//	m->z[i] = m->k[i] * m->r[i];

		dot_kernel<<<grid, block>>>(dot_d, z_d, z_d, n);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
		pnorm = sqrt(dot_h);
		//pnorm = sqrt(get_dot(m->z, m->z, m->nrow));

		dot_kernel<<<grid, block>>>(dot_d, z_d, z_d, n);
		cudaMemcpy(&dot_h, dot_d, sizeof(double), cudaMemcpyDeviceToHost);
		double rz_n = dot_h;
		//double rz_n = get_dot(m->r, m->z, m->nrow);

		const double beta = rz_n / rz;
	        axpby<<<grid, block>>>(1.0, z_d, beta, p_d, n);
		//for (int i = 0; i < m->nrow; ++i)
		//	m->p[i] = m->z[i] + beta * m->p[i];

		rz = rz_n;
		its++;
	}

	*err = rz;

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
