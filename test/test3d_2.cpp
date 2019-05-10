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
#include <cassert>
#include <chrono>
#include <bits/stdc++.h>


#include "micro.hpp"


using namespace std;
using namespace std::chrono;


#define D_EPS 5.0e-4


int main(int argc, char **argv)
{

	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " n [ngp] [steps]" << endl;
		return(1);
	}

	const int dir = 1;
	const int n = atoi(argv[1]);
	const int ngp = (argc > 2 ? atoi(argv[2]) : 10);
	const int time_steps = (argc > 3 ? atoi(argv[3]) : 10);  // Optional value

	assert(n > 1 && ngp > 0 && time_steps > 0);

	const int size[3] = { n, n, n };
	const int micro_type = MIC_SPHERE; // 2 materiales matriz y fibra (3D esfera en matriz)
	const double micro_params[4] = { 1.0, 1.0, 1.0, 0.1 };

	material_base mat_params[2];
	material_set(&mat_params[0], 0, 1.0e7, 0.3, 0.0, 0.0, 1.0e1);
	material_set(&mat_params[1], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);

	micropp<3> micro(ngp, size, micro_type, micro_params, mat_params, ONE_WAY, true, 5);
	micro.print_info();

	auto start = high_resolution_clock::now();

	double sig[6];
	double eps[6] = { 0. };

	cout << scientific;

	for (int t = 0; t < time_steps; ++t) {

		cout << "Time step = " << t << endl;

		if (t < 20)
			eps[dir] += D_EPS;
		else if (t < 40)
			eps[dir] -= D_EPS;
		else
			eps[dir] += D_EPS;

		for (int gp = 0; gp < ngp; ++gp) {
			micro.set_strain(gp, eps);
		}
		cout << "eps = ";
		for (int i = 0; i < 6; ++i) {
			cout << eps[i] << " ";
		}
		cout << endl;

		cout << "Homogenizing ..." << endl;
		micro.homogenize();

		for (int gp = 0; gp < ngp; ++gp) {
			cout << "sig = ";
			micro.get_stress(gp, sig);
			for (int i = 0; i < 6; ++i) {
				cout << sig[i] << "\t";
			}
			cout << endl;
		}
		cout << endl;

		micro.update_vars();
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "time = " << duration.count() << " ms" << endl;

	return 0;
}
