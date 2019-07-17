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

using namespace std;


double eps_vs_t(double time, double t_final)
{
	double eps;
	const double eps_max = 1.0e-1;
	const double t1 = (1.0 / 4.0) * t_final;
	const double t2 = (2.0 / 4.0) * t_final;
	const double t3 = (3.0 / 4.0) * t_final;
	if (time < t1) {
		eps = 0.5 * eps_max * time;
	} else if (time >= t1 && time < t2) {
		eps = 0.5 * eps_max * (t2 - time);
	} else if (time >= t2 && time < t3) {
		eps = eps_max * (time - t2);
	} else {
		eps = eps_max * (t_final - time);
	}
	return eps;
}


int main (int argc, char *argv[])
{
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " n [steps]" << endl;
		return(1);
	}

	const int n = atoi(argv[1]);

	const int print = (argc > 2) ? atoi(argv[2]) : 0;
	if (print < 0 || print > 1) {
		cerr << "Error in [print] argument, it only can be 0 or 1" << endl;
		return(1);
	}

	const int time_steps = (argc > 3 ? atoi(argv[3]) : 10);

	const int dir = 1;
	const double t_final = 1.0;
	const double dt = t_final / time_steps;
	double time = 0.0;

	ofstream file;
	file.open("result.dat");

	micropp_params_t mic_params;

	mic_params.size[0] = n;
	mic_params.size[1] = n;
	mic_params.size[2] = n;
	mic_params.type = MIC_HOMOGENEOUS;
	mic_params.nr_max_its = 8;
	material_set(&mic_params.materials[0], 2, 1.0e7, 0.3, 0.0, 0.0, 1.0e5);
	material_set(&mic_params.materials[1], 0, 0.0, 0.0, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[2], 0, 0.0, 0.0, 0.0, 0.0, 0.0);
	mic_params.lin_stress = false;

	mic_params.print();

	micropp<3> micro(mic_params);
	micro.print_info();

	double sig[6];
	double ctan[36];
	double eps[6] = { 0. };

	cout << scientific;

	for (int t = 0; t < time_steps; ++t) {

		cout << "time step = " << t << endl;

		eps[dir] = eps_vs_t(time, t_final);

		micro.set_strain(0, eps);
		micro.homogenize();
		micro.get_stress(0, sig);
		micro.get_ctan(0, ctan);
		int non_linear = micro.is_non_linear(0);
		int cost = micro.get_cost(0);
		bool has_converged = micro.has_converged(0);

		char filename[128];
		snprintf(filename, 128, "test_damage_%d", t);

		if (print) {
			char filename[128];
			snprintf(filename, 128, "micropp_%d", t);
			micro.output (0, filename);
		}

		micro.update_vars();

		cout 	<< "NL        = " << non_linear << endl
			<< "Cost      = " << cost << endl
			<< "Converged = " << has_converged << endl;

		cout << "eps =\t";
		for (int i = 0; i < 6; ++i)
			cout << eps[i] << "\t";
		cout << endl;

		cout << "sig =\t";
		for (int i = 0; i < 6; ++i)
			cout << sig[i] << "\t";
		cout << endl;

		file    << setw(14)
			<< eps[dir] << "\t"
			<< sig[dir] << "\t" << endl;

		time += dt;
	}

	file.close();
	return 0;
}
