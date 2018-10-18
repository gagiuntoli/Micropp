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

#include "micropp_c_wrapper.h"

#define D_EPS 5.0e-4

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
	const int time_steps = (argc > 5 ? atoi(argv[5]) : 10);  // Optional value
	int size[3] = { nx, ny, nz };

    micropp_C_material_create();
    micropp_C_material_set(0, 1.0e7, 0.25, 1.0e4, 1.0e7, 0);
    micropp_C_material_set(1, 1.0e7, 0.25, 1.0e4, 1.0e7, 0);
    micropp_C_material_print(0);
    micropp_C_material_print(1);

    int ngpl = 1;
    int type = 1;
    double params[4] = { 1., 1., 1., .5 };
    micropp_C_create3(ngpl, size, type, params);
    micropp_C_print_info();


	double sig[6], ctan[36];
	double eps[6] = { 0. };
    int i, j, t;

	for (t = 0; t < time_steps; ++t) {

        printf("time step = %d\n", t);

		if (t < 80)
			eps[dir] += D_EPS;
		else if (t < 160)
			eps[dir] -= D_EPS;
		else
			eps[dir] += D_EPS;

		micropp_C_set_strain3(0, eps);
		micropp_C_homogenize();
		micropp_C_get_stress3(0, sig);
		micropp_C_get_ctan3(0, ctan);

		micropp_C_update_vars();

        printf("\neps =\n");
		for (i = 0; i < 6; ++i)
			printf("%e\t", eps[i]);
        printf("\n");

        printf("\nsig =\n");
		for (i = 0; i < 6; ++i)
			printf("%e\t", sig[i]);
        printf("\n");

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
