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
#include <chrono>
#include <bits/stdc++.h>

#include "micro.hpp"

using namespace std;
using namespace std::chrono;

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

	int dir = 2;
	double eps[nvoi] = { 0.0 };
	double sig[nvoi];

	micropp_params_t mic_params;

	mic_params.ngp = ngp;
	mic_params.size[0] = n;
	mic_params.size[1] = n;
	mic_params.size[2] = n;
	mic_params.type = MIC_SPHERE;

	material_set(&mic_params.materials[0], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[1], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[2], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	mic_params.calc_ctan_lin = false;
	mic_params.use_A0 = false;
	mic_params.lin_stress = false;

	mic_params.print();

	micropp<3> micro(mic_params);
	micro.print_info();

	auto start = high_resolution_clock::now();

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
			micro.set_strain(gp, eps);
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
			micro.get_stress(gp, sig);
			cout << "gp = " << gp << " sig  = ";
			for (int i = 0; i < 6; ++i)
				cout << setw(14) << sig[i] << "\t";
			cout << endl;
		}

		micro.update_vars();

	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "time = " << duration.count() << " ms" << endl;

	return 0;
}
