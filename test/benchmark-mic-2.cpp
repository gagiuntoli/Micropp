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


#include "micropp.hpp"


using namespace std;
using namespace std::chrono;


double eps_vs_t(double time, double t_final) {
	const double eps_max = 1.0e-1;
	double eps = eps_max * time;
	return eps;
}


int main(int argc, char **argv)
{

	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " n [print=0|1] [steps]" << endl;
		return(1);
	}

	const int n = atoi(argv[1]);
	const int ngp = 1;
	const int print = (argc > 2 ? atoi(argv[2]) : 0);  // Optional value
	const int time_steps = (argc > 3 ? atoi(argv[3]) : 10);  // Optional value

	const int dir = 0;
	const double t_final = 0.15;
	const double dt = t_final / time_steps;
	double time = 0.0;

	assert(n > 1 && 0 <= print && print < 2 && time_steps >= 0);

	micropp_params_t mic_params;

	mic_params.ngp = ngp;
	mic_params.size[0] = n;
	mic_params.size[1] = n;
	mic_params.size[2] = n;
	mic_params.type = MIC3D_SPHERES;
	mic_params.geo_params[0] = 0.0;
	mic_params.geo_params[1] = 0.0;
	mic_params.geo_params[2] = 0.0;
	mic_params.geo_params[3] = 0.0;
	//mic_params.subiterations = true;
	//mic_params.nsubiterations = 10;
	mic_params.nr_max_its = 12;
	material_set(&mic_params.materials[0], 0, 1.0e7, 0.25, 1.0e7, 1.0e7, 0.0);
	material_set(&mic_params.materials[1], 0, 2.0e7, 0.25,   0.0,   0.0, 0.0);
	material_set(&mic_params.materials[2], 0, 3.0e7, 0.3,   0.0,   0.0, 0.0);
	mic_params.mpi_rank = 0;
	mic_params.calc_ctan_lin = false;
	mic_params.lin_stress = false;

	mic_params.print();

	micropp<3> micro(mic_params);
	//micro.print_info();

	if (print) {
		char filename[128];
		snprintf(filename, 128, "micropp_%d", 0);
		micro.output (0, filename);
	}

	ofstream file;
	file.open("result.dat");

	auto start = high_resolution_clock::now();

	double sig[6];
	double eps[6] = { 0. };

	cout << scientific;

	for (int t = 0; t < time_steps; ++t) {

		cout << "Time step = " << t << endl;

		eps[dir] = eps_vs_t(time, t_final);

		micro.set_strain(0, eps);

		cout << "eps = ";
		for (int i = 0; i < 6; ++i) {
			cout << eps[i] << " ";
		}
		cout << endl;

		cout << "Homogenizing ..." << endl;
		micro.homogenize();

		int non_linear = micro.is_non_linear(0);
		int cost = micro.get_cost(0);
		bool has_converged = micro.has_converged(0);

		cout 	<< "NL        = " << non_linear << endl
			<< "Cost      = " << cost << endl
			<< "Converged = " << has_converged << endl;

		cout << "sig = ";
		micro.get_stress(0, sig);
		for (int i = 0; i < 6; ++i) {
			cout << sig[i] << "\t";
		}
		cout << endl;

		micro.update_vars();

		file    << setw(14)
			<< eps[dir] << "\t"
			<< sig[dir] << "\t" << endl;

		if (print) {
			char filename[128];
			snprintf(filename, 128, "micropp_%d", t);
			micro.output (0, filename);
		}

		time += dt;
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "time = " << duration.count() << " ms" << endl;

	return 0;
}
