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


#include <stdio.h>
#include <stdlib.h>


#include "micropp_c.h"

#define D_EPS 2.0e-4


int main (int argc, char *argv[])
{
	const int dim = 3;
	if (argc < 5) {
		printf("Usage: %s nx ny nz dir [steps]\n", argv[0]);
		return(1);
	}

	const int nx = atoi(argv[1]);
	const int ny = atoi(argv[2]);
	const int nz = atoi(argv[3]);
	const int dir = atoi(argv[4]);
	const int time_steps = (argc > 5 ? atoi(argv[5]) : 10); // Optional
	int size[3] = { nx, ny, nz };

	struct material_base matlist[2];
	material_set(&matlist[0], 1, 1.0e7, 0.25, 1.0e4, 1.0e4, 1.0);
	material_set(&matlist[1], 0, 1.0e7, 0.25, 1.0e4, 1.0e7, 0.0);
	material_print(&matlist[0]);
	material_print(&matlist[1]);

	struct micropp3 micro;
	int ngpl = 1;
	int micro_type = 5;
	int nsubiterations = 1;
	int mpi_rank = 0;
	double params[4] = { 1., 1., 1., .15 };
	micropp3_new(&micro, ngpl, size, micro_type, params, matlist, nsubiterations, mpi_rank);
	micropp3_print_info(&micro);

	double sig[6], ctan[36];
	double eps[6] = { 0. };
	int i, j, t;

	for (t = 0; t < time_steps; ++t) {

		printf("\nTime step = %d\n", t);

		if (t < 80)
			eps[dir] += D_EPS;
		else if (t < 160)
			eps[dir] -= D_EPS;
		else
			eps[dir] += D_EPS;

		micropp3_set_strain(&micro, 0, eps);
		micropp3_homogenize(&micro);
		micropp3_get_stress(&micro, 0, sig);
		micropp3_get_ctan(&micro, 0, ctan);
		int sigma_cost = micropp3_get_cost(&micro, 0);
		bool non_linear = micropp3_is_non_linear(&micro, 0);
		int num_non_linear = micropp3_get_non_linear_gps(&micro);

		micropp3_update_vars(&micro);

		char filename[128];
		sprintf(filename, "test3d_6_%d", t);
		micropp3_output(&micro, 0, filename);

		printf("sigma_cost       = %d\n", sigma_cost);
		printf("Non-Linear       = %d\n", non_linear);
		printf("Non-Linear Total = %d\n", num_non_linear);

		printf("eps =\n");
		for (i = 0; i < 6; ++i)
			printf("%e\t", eps[i]);

		printf("\nsig =\n");
		for (i = 0; i < 6; ++i)
			printf("%e\t", sig[i]);

		printf("\nctan =\n");
		for (i = 0; i < 6; ++i) {
			for (j = 0; j < 6; ++j)
				printf("%e\t", ctan[i * 6 + j]);
			printf("\n");
		}
		printf("\n");

	}
	return 0;
}
