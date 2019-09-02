/*
 *  This is a test example for MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
 *                         Judicaël Grasset <judicael.grasset@stfc.ac.uk>
 *                         Alejandro Figueroa <afiguer7@maisonlive.gmu.edu>
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

#include "micropp.hpp"

using namespace std;

#define D_EPS 5.0e-4

int main (int argc, char *argv[])
{
	// Execution ./test3d_1 n [print] [steps]
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " n [use_A0 = 0|1] [steps]" << endl;
		return(1);
	}

	const int dir = 1;
	const int n = atoi(argv[1]);

	const int use_A0 = (argc > 2) ? atoi(argv[2]) : 0;
	if (use_A0 < 0 || use_A0 > 1) {
		cerr    << "Error in [use_A0] argument, it only can be 0 or 1"
			<< endl;
		exit(1);
	}

	const int time_steps = (argc > 3 ? atoi(argv[3]) : 10);

	micropp_params_t mic_params;

	mic_params.ngp = 1;
	mic_params.size[0] = n;
	mic_params.size[1] = n;
	mic_params.size[2] = n;
	mic_params.type = MIC_SPHERE;
	material_set(&mic_params.materials[0], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[1], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[2], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	mic_params.calc_ctan_lin = true;
	mic_params.use_A0 = use_A0;

	mic_params.print();

	micropp<3> micro(mic_params);

	micro.print_info();

	double sig[6];
	double ctan[36];
	double eps[6] = { 0. };

	cout << scientific;

	for (int t = 0; t < time_steps; ++t) {

		cout << "time step = " << t << endl;

		if (t < 80)
			eps[dir] += D_EPS;
		else if (t < 160)
			eps[dir] -= D_EPS;
		else
			eps[dir] += D_EPS;

		micro.set_strain(0, eps);
		micro.homogenize();
		micro.get_stress(0, sig);
		micro.get_ctan(0, ctan);
		int non_linear = micro.is_non_linear(0);
		int cost = micro.get_cost(0);
		bool has_converged = micro.has_converged(0);

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

		cout << "ctan = " << endl;
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j) {
				cout << ctan[i * 6 + j] << "\t";
			}
			cout << endl;
		}
		cout << endl;
	}

	return 0;
}
