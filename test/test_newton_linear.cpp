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

#define REPETITIONS 5

using namespace std;

const double strain[6] = { 1., 2., 3., 1., 1., 1. };
double s;

int main (int argc, char *argv[])
{
	class test_t : public micropp<3> {

		public:
			test_t(const int size[3], const int micro_type, const double micro_params[5],
			       const material_t mat_params[2])
				:micropp<3> (1, size, micro_type, micro_params, mat_params, &s)
			{};

			~test_t() {};

			void just_do_it(void)
			{

#pragma omp parallel for
				for (int i = 0; i < REPETITIONS; ++i) {
					int thread_id = omp_get_thread_num();
					memset(u[thread_id], 0.0, nndim * sizeof(double));
					newton_raphson_linear(&A0[thread_id],
							      b[thread_id],
							      u[thread_id],
							      du[thread_id],
							      strain,
							      true);
					cout << endl;
				}

			};

	};

	if (argc < 2) {
		/* argv[1] (n) : Problem size
		 * argv[2] (a) : Factor that says Ef = Em x a
		 */
		cerr << "Usage: " << argv[0] << " [n = 10] [a = 1]" << endl;
	}

	const int n = (argc > 1) ? atoi(argv[1]) : 10;
	const double a = (argc > 2) ? atoi(argv[2]) : 1.0;

	int size[3] = { n, n, n };

	int micro_type = 2;
	double micro_params[4] = { 1., 1., 1., 0.2 };

	double Em = 1.0e8;

	material_t mat_params[2];
	mat_params[0].set(Em, 0.25, 1.0e8, 1.0e4, 0);
	mat_params[1].set(Em * a, 0.25, 1.0e8, 1.0e4, 0);

	test_t test(size, micro_type, micro_params, mat_params);
	test.just_do_it();

	return 0;
}
