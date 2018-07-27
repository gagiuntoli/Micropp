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

#include "ell.hpp"

using namespace std;

int main (int argc, char *argv[])
{
	if (argc < 3) {
		cerr << "Usage: " << argv[0] << " nx ny [nz]" << endl;
		return(1);
	}

	const int nx = atoi(argv[1]);
	const int ny = atoi(argv[2]);
	const int nz = (argc > 3) ? atoi(argv[3]) : 1;
	const int dim = (argc > 3) ? 3 : 2;

	assert(nx > 1 && ny > 1 && nz > 0);

	const int nex = nx - 1;
	const int ney = ny - 1;
	const int nez = (dim == 3) ? nz - 1 : 0;
	const int ns[3] = { nx, ny, nz };
	const int nfield = 1;

	ell_matrix A1;

	int cg_its;
	double cg_err;

	if (dim == 2) {

		ell_init(&A1, nfield, dim, ns, 1.0e-5, 20);

		cout << "A1.nrow =\t" << A1.nrow << endl;
		cout << "A1.ncol =\t" << A1.ncol << endl;
		cout << "A1.nnz =\t" << A1.nnz << endl;
		assert( A1.nrow == (nx * ny) && 
				A1.ncol == (nx * ny) && 
				A1.nnz == nfield * 9 );

		const double Ae[4 * 4] = {
			10.6667,   -2.6667,   -5.3333,   -2.6667,
			-2.6667,   10.6667,   -2.6667,   -5.3333,
			-5.3333,   -2.6667,   10.6667,   -2.6667,
			-2.6667,   -5.3333,   -2.6667,   10.6667 };

		ell_set_zero_mat(&A1);
		for (int ex = 0; ex < nex; ++ex)
			for (int ey = 0; ey < ney; ++ey)
				ell_add_2D(&A1, ex, ey, Ae);

		ell_set_bc_2D(&A1);
		for (int i = 0; i < A1.nrow; ++i) {
			for (int j = 0; j < 9; ++j)
				cout << A1.vals[i * A1.nnz + j] << " ";
			cout << endl;
		}

		double *x = (double *)malloc(A1.nrow * sizeof(double));
		double *b = (double *)calloc(A1.nrow, sizeof(double));
		for (int i = 0; i < A1.nrow; ++i)
			b[i] = 1.0;

		ell_solve_cgpd(&A1, b, x, &cg_err, &cg_its);
		cout << "Err =\t" << cg_err << "\tIts =\t" << cg_its << endl;

	} else if (dim == 3) {

		ell_init(&A1, nfield, dim, ns, 1.0e-5, 20);

		cout << "A1.nrow =\t" << A1.nrow << endl;
		cout << "A1.ncol =\t" << A1.ncol << endl;
		cout << "A1.nnz =\t" << A1.nnz << endl;
		assert( A1.nrow == (nx * ny * nz) && 
				A1.ncol == (nx * ny * nz) && 
				A1.nnz == nfield * 27 );

	}

	ell_free(&A1);

	return 0;
}
