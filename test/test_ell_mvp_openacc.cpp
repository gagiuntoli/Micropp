/*
 *  This is a test example for MicroPP: a finite element library
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

#include <iostream>
#include <iomanip>

#include <ctime>
#include <cassert>
#include <chrono>
#include <cmath>

#include "ell.hpp"

using namespace std;
using namespace std::chrono; 


void mvp_gpu_1(const int *cols, const double *vals, const int nnz, double *y, const double *x, const int nn)
{
#pragma acc parallel loop gang worker present (cols[:nn * nnz], vals[:nn * nnz], x[:nn], y[:nn])
	for (int i = 0; i < nn; i++) {
		double tmp = 0.0;
		const int ix = i * nnz;
#pragma acc loop vector
		for (int j = 0; j < nnz; j++) {
			const int index = ix + j;
			const int col_id = cols[index];
			tmp += vals[index] * x[col_id];
		}
		y[i] = tmp;
	}
}


void mvp_gpu_2(const ell_matrix *A, double *y, const double *x)
{
	const int nn = A->nrow;
	const int nnz = A->nnz;
	const int *cols = A->cols;
	const double *vals = A->vals;
#pragma acc parallel loop gang worker present (cols[:nn * nnz], vals[:nn * nnz], x[:nn], y[:nn])
	for (int i = 0; i < nn; i++) {
		double tmp = 0.0;
		const int ix = i * nnz;
#pragma acc loop vector
		for (int j = 0; j < nnz; j++) {
			const int index = ix + j;
			const int col_id = cols[index];
			tmp += vals[index] * x[col_id];
		}
		y[i] = tmp;
	}
}


void mvp_gpu_3(const ell_matrix *A, double *y, const double *x)
{
	const int nn = A->nrow;
	const int nnz = A->nnz;
#pragma acc parallel loop gang worker present (A[:1], A->cols[:nn * nnz], A->vals[:nn * nnz], x[:nn], y[:nn])
	for (int i = 0; i < nn; i++) {
		double tmp = 0.0;
		const int ix = i * nnz;
#pragma acc loop vector
		for (int j = 0; j < nnz; j++) {
			const int index = ix + j;
			const int col_id = A->cols[index];
			tmp += A->vals[index] * x[col_id];
		}
		y[i] = tmp;
	}
}


int main (int argc, char *argv[])
{
	const int n = (argc > 1) ? atoi(argv[1]) : 10;
	const int nmul = (argc > 2) ? atoi(argv[2]) : 10;
	const int dim = 3;
	const int size[3] = { n + 1 , n + 1, n + 1 };
	const int nfield = 3;

	assert(n > 0 && nmul > 0);

	ell_matrix A1;
	ell_init(&A1, nfield, dim, size, CGPD, 1.0e-5, 1.0e-5, 20);

	cout << "Matrix size = " << A1.nrow << " x " << A1.nrow << endl;
	cout << "Number of non-zeros per row = " << A1.nnz << endl;
	cout << "Number of multiplications = " << nmul << endl;
	const int nn = A1.nrow;

	double Ae[24 * 24];
	for (int i = 0; i < 24; ++i)
		for (int j = 0; j < 24; ++j)
			Ae[i * 24 + j] = i + j;

	ell_set_zero_mat(&A1);
	for (int ex = 0; ex < n; ++ex)
		for (int ey = 0; ey < n; ++ey)
			for (int ez = 0; ez < n; ++ez)
				ell_add_3D(&A1, ex, ey, ez, Ae);

	ell_set_bc_3D(&A1);

	double *x = (double *)malloc(nn * sizeof(double));
	double *y = (double *)malloc(nn * sizeof(double));

	double normx = 0.0;
	double normy = 0.0;
	int nnz = A1.nnz;

	for (int i = 0; i < nn; ++i)
		x[i] = 1.0;

#pragma acc data copyin(x[:nn])
	{
		auto start = high_resolution_clock::now();

		cout << "Test mvp_gpu_1" << endl;
		const int *cols = A1.cols;
		const double *vals = A1.vals;

#pragma acc data copyin(cols[:nnz * nn], vals[:nnz * nn])
#pragma acc data copy(y[:nn])
		for (int its = 0; its < nmul; its++) {

			mvp_gpu_1(cols, vals, nnz, y, x, nn);

			normy = 0.0;
#pragma acc parallel loop
			for (int i = 0; i < nn; i++) {
				normy += y[i] * y[i];
			}
			normy = sqrt(normy);
		}
		cout << "|x| = " << normx << " |y| = " << normy << endl;
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		cout << "mpv total time = " << duration.count() << " ms" << endl; 


		start = high_resolution_clock::now();
		cout << endl << "Test mvp_gpu_2" << endl;

#pragma acc data copy(x[:nn])
#pragma acc data copyin(A1)
#pragma acc data copyin(A1.cols[:nnz * nn], A1.vals[:nnz * nn])
#pragma acc data copy(y[:nn])
		for (int its = 0; its < nmul; its++) {

			mvp_gpu_2(&A1, y, x);

			normy = 0.0;
#pragma acc parallel loop
			for (int i = 0; i < nn; i++) {
				normy += y[i] * y[i];
			}
			normy = sqrt(normy);
		}
		cout << "|x| = " << normx << " |y| = " << normy << endl;
		stop = high_resolution_clock::now();
		duration = duration_cast<milliseconds>(stop - start);
		cout << "mpv total time = " << duration.count() << " ms" << endl; 


		start = high_resolution_clock::now();
		cout << endl << "Test mvp_gpu_3" << endl;

#pragma acc data copy(x[:nn])
#pragma acc data copyin(A1)
#pragma acc data copyin(A1.cols[:nnz * nn], A1.vals[:nnz * nn])
#pragma acc data copy(y[:nn])
		for (int its = 0; its < nmul; its++) {

			mvp_gpu_3(&A1, y, x);

			normy = 0.0;
#pragma acc parallel loop
			for (int i = 0; i < nn; i++) {
				normy += y[i] * y[i];
			}
			normy = sqrt(normy);
		}
		cout << "|x| = " << normx << " |y| = " << normy << endl;
		stop = high_resolution_clock::now();
		duration = duration_cast<milliseconds>(stop - start);
		cout << "mpv total time = " << duration.count() << " ms" << endl; 
	}

	ell_free(&A1);
	free(x);
	free(y);

	return 0;
}
