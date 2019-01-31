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

#define REPETITIONS 20

using namespace std;


int main (int argc, char *argv[])
{
	class test_t : public micropp<3> {

		public:
			test_t(const int size[3], const int micro_type, const double micro_params[5],
			       const material_t mat_params[2])
				:micropp<3> (1, size, micro_type, micro_params, mat_params, NO_COUPLING)
			{};

			~test_t() {};

			void just_do_it(const int size[3])
			{
				int nthreads = omp_get_max_threads();
				double **u = (double**)malloc(nthreads * sizeof(double*));
				ell_matrix *A = (ell_matrix*)malloc(nthreads * sizeof(ell_matrix));
				for (int i = 0; i < nthreads; ++i) {
					u[i] = (double*)malloc(nndim * sizeof(double));
					ell_init(&A[i], 3, 3, size, CG_MIN_ERR, CG_REL_ERR, CG_MAX_ITS);
				}

#pragma omp parallel for
				for (int i = 0; i < REPETITIONS; ++i) {
					const double strain[6] = { 1., 2., 3., 1., 1., 1. };
					int thread_id = omp_get_thread_num();
					memset(u[thread_id], 0.0, nndim * sizeof(double));
					set_displ_bc(strain, u[thread_id]);
					assembly_mat(&A[thread_id], u[thread_id], nullptr);
					double norm = ell_get_norm(&A[thread_id]);
					printf("Thread id : %d norm = %lf\n", thread_id, norm);
				}

				for (int i = 0; i < nthreads; ++i) {
					free(u[i]);
				}
				free(u);
			};

	};

	if (argc < 2) {
		/* argv[1] (n) : Problem size
		 */
		cerr << "Usage: " << argv[0] << " [n = 10]" << endl;
	}

	const int n = (argc > 1) ? atoi(argv[1]) : 10;

	int size[3] = { n, n, n };

	int micro_type = 2;
	double micro_params[4] = { 1., 1., 1., 0.2 };

	double Em = 1.0e8;

	material_t mat_params[2];
	mat_params[0].set(Em, 0.25, 1.0e8, 1.0e4, 0);
	mat_params[1].set(Em, 0.25, 1.0e8, 1.0e4, 0);

	printf("Initializing\n");
	test_t test(size, micro_type, micro_params, mat_params);

	printf("Doing test\n");
	double time = omp_get_wtime();
	test.just_do_it(size);
	time = omp_get_wtime() - time;
	printf("time = %lf\n", time);

	return 0;
}
