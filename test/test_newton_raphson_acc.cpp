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


#include "micro.hpp"


#define REPETITIONS 1


using namespace std;


class test_t : public micropp<3> {

	public:
		test_t(const int size[3], const int micro_type, const double micro_params[5],
		       const material_base materials[2])
			:micropp<3> (1, size, micro_type, micro_params, materials, NO_COUPLING)
		{};

		~test_t() {};

		void newton_raphson(const double *strain)
		{
			const int ns[3] = { nx, ny, nz };
			ell_matrix A;
			ell_init(&A, dim, dim, ns);

			double *b = (double *) calloc(nndim, sizeof(double));
			double *du = (double *) calloc(nndim, sizeof(double));
			double *u = (double *) calloc(nndim, sizeof(double));

			newton_t newton = newton_raphson_acc(&A, b, u, du, strain, nullptr);
			newton.print();

			ell_free(&A);
			free(b);
			free(u);
			free(du);

		};

};

int main (int argc, char *argv[])
{
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " [n = 10] " << endl;
		exit(1);
	}

	const int n = (argc > 1) ? atoi(argv[1]) : 10;
	const int size[3] = { n, n, n };
	const int micro_type = 2;
	const double micro_params[4] = { 1., 1., 1., 0.2 };

	material_base mat_params[2];
	material_set(&mat_params[0], 0, 1.0e7, 0.3, 0.0, 0.0, 1.0e1);
	material_set(&mat_params[1], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);

	const double strain[6] = { 1., 2., 3., 1., 1., 1. };

	test_t test(size, micro_type, micro_params, mat_params);

	for (int i = 0; i < REPETITIONS; ++i)
		test.newton_raphson(strain);

	return 0;
}
