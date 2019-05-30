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
		test_t(const micropp_params_t params)
			:micropp<3> (params)
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

			newton_t newton = newton_raphson_acc(&A, b, u, du, strain);
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
	const double strain[6] = { 1., 2., 3., 1., 1., 1. };

	micropp_params_t mic_params;

	mic_params.ngp = 1;
	mic_params.size[0] = n;
	mic_params.size[1] = n;
	mic_params.size[2] = n;
	mic_params.type = MIC_SPHERE;
	material_set(&mic_params.materials[0], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[1], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[2], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);

	mic_params.print();

	test_t test(mic_params);

	for (int i = 0; i < REPETITIONS; ++i)
		test.newton_raphson(strain);

	return 0;
}
