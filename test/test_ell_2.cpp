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
	const int nx = 3;
	const int ny = 3;
	const int nz = 3;

	ell_matrix A1;

	const int ns[3] = { nx, ny, nz };
	const int nfield = 1;
	const int dim = 3;
	ell_init(&A1, nfield, dim, ns, 1.0e-5, 1.0e-5, 20);

	cout << "A1.nrow =\t" << A1.nrow << endl;
	cout << "A1.ncol =\t" << A1.ncol << endl;
	cout << "A1.nnz =\t" << A1.nnz << endl;
	assert( A1.nrow == (nx * ny * nz) &&
		A1.ncol == (nx * ny * nz) &&
		A1.nnz == nfield * 27 );

	return 0;
}
