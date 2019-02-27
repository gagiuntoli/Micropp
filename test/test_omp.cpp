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

#include <cstring>
#include <ctime>
#include <cassert>

#include <bits/stdc++.h>

#include "micro.hpp"

using namespace std;

#define D_EPS 0.01

int main(int argc, char **argv)
{

	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " n [ngp] [steps]" << endl;
		return(1);
	}

	const int nvoi = 6;
	const int n = atoi(argv[1]);
	const int ngp = (argc > 2 ? atoi(argv[2]) : 2);
	const int time_steps = (argc > 3 ? atoi(argv[3]) : 1);  // Optional value

	assert(n > 1 && ngp > 0 && time_steps > 0);

	int size[3] = { n, n, n };

	const int micro_type = 3;
	const double micro_params[4] = { 1.0, 1.0, 1.0, 0.1 };

	material_t mat_params[2];
	mat_params[0].set(1.0e6, 0.3, 5.0e4, 5.0e4, 0);
	mat_params[1].set(1.0e6, 0.3, 1.0e4, 0.0e0, 0);

	int dir = 2;
	double eps[nvoi] = { 0.0 };
	double sig[nvoi];
	double ctan[nvoi * nvoi];

	micropp<3> micro(ngp, size, micro_type, micro_params, mat_params, NO_COUPLING);

	double time;
#ifdef _OPENMP
	time = omp_get_wtime() - time;
#else
	time = clock() - time;
#endif

	cout << scientific;
	for (int t = 0; t < time_steps; ++t) {

		cout << "Time step = " << t << endl;

		if (t < 20)
			eps[dir] += D_EPS;
		else if (t < 40)
			eps[dir] -= D_EPS;
		else
			eps[dir] += D_EPS;

		cout << "setting strains ..." << endl;
		for (int gp = 0; gp < ngp; ++gp) {
			micro.set_macro_strain(gp, eps);
			cout << "gp = " << gp << " eps = ";
			cout << scientific;

			for (int i = 0; i < 6; ++i)
				cout << setw(14) << eps[i] << " ";
			cout << endl;
		}

		cout << "Homogenizing ..." << endl;
		micro.homogenize(); // This function uses OMP

		cout << "Getting stresses ..." << endl;
		for (int gp = 0; gp < ngp; ++gp) {
			micro.get_macro_stress(gp, sig);
			cout << "gp = " << gp << " sig  = ";
			for (int i = 0; i < 6; ++i)
				cout << setw(14) << sig[i] << "\t";
			cout << endl;
		}

		micro.update_vars();

	}

#ifdef _OPENMP
	time = omp_get_wtime() - time;
#else
	time = clock() - time;
#endif
	printf("time = %lf\n", time);

	return 0;
}
