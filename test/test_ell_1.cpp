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
#include <cmath>
#include <cassert>

#include "ell.hpp"

using namespace std;

int main (int argc, char *argv[])
{
	const int nx = 3;
	const int ny = 3;
	const int nex = nx - 1;
	const int ney = ny - 1;

	ell_matrix A1;
	const double Ae[4 * 4] = {
		10.6667,   -2.6667,   -5.3333,   -2.6667,
		-2.6667,   10.6667,   -2.6667,   -5.3333,
		-5.3333,   -2.6667,   10.6667,   -2.6667,
		-2.6667,   -5.3333,   -2.6667,   10.6667 };

	double cg_err;

	const int ns[3] = { nx, ny, 0 };
	const int nfield = 1;
	const int dim = 2;
	ell_init(&A1, nfield, dim, ns, 1.0e-5, 50);

	cout << "A1.nrow =\t" << A1.nrow << endl;
	cout << "A1.ncol =\t" << A1.ncol << endl;
	cout << "A1.nnz =\t" << A1.nnz << endl;
	assert( A1.nrow == (nx * ny) && 
			A1.ncol == (nx * ny) && 
			A1.nnz == nfield * 9 );

    const int cols_sol[9 * 9] = {
		0, 0, 0, 0, 0, 1, 0, 3, 4,
		0, 0, 0, 0, 1, 2, 3, 4, 5,
		0, 0, 0, 1, 2, 0, 4, 5, 0,
		0, 0, 1, 0, 3, 4, 0, 6, 7,
		0, 1, 2, 3, 4, 5, 6, 7, 8,
		1, 2, 0, 4, 5, 0, 7, 8, 0,
		0, 3, 4, 0, 6, 7, 0, 0, 0,
		3, 4, 5, 6, 7, 8, 0, 0, 0,
		4, 5, 0, 7, 8, 0, 0, 0, 0 };
   
	for (int i = 0; i < 9 * 9; ++i)
		assert (A1.cols[i] == cols_sol[i]);

	assert( A1.nrow == (nx * ny) && 
			A1.ncol == (nx * ny) && 
			A1.nnz == nfield * 9 );

	ell_set_zero_mat(&A1);
	for (int ex = 0; ex < nex; ++ex)
		for (int ey = 0; ey < ney; ++ey)
			ell_add_2D(&A1, ex, ey, Ae);

	ell_set_bc_2D(&A1);

	const double dir_nod_sol[9] = { 0, 0, 0, 0, 1, 0, 0, 0, 0 };
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 9; ++j)
			assert (A1.vals[i * A1.nnz + j] == dir_nod_sol[j]);

	for (int i = 5; i < 9; ++i)
		for (int j = 0; j < 9; ++j)
			assert (A1.vals[i * A1.nnz + j] == dir_nod_sol[j]);

	const double int_nod_sol[9] = {
		-5.3333, -5.3333, -5.3333, -5.3333, 42.6668, -5.3333, -5.3333, -5.3333, -5.3333 };
	for (int j = 0; j < 9; ++j)
		assert (fabs(A1.vals[4 * A1.nnz + j] - int_nod_sol[j]) < 1.0e-2);

	double *x = (double *)calloc(A1.nrow, sizeof(double));
	double *b = (double *)calloc(A1.nrow, sizeof(double));
	b[4] = 1.0;

	int cg_its = ell_solve_cgpd(&A1, b, x, &cg_err);

	const double x_sol[9] = { 0, 0, 0, 0, 0.0234374, 0, 0, 0, 0 };
	for (int j = 0; j < 9; ++j)
		assert (fabs(x[j] - x_sol[j]) < 1.0e-3);

	assert (cg_err < 1.0e-10);
	assert (cg_its == 1);

	ell_free(&A1);

	return 0;
}
