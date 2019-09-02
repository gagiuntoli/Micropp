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

const int time_steps = 10;
const double sig_sol[time_steps][6] = {
	{ 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.0, 0.0, 0.0 },
	{ 5.40000000000000e+04, 1.80000000000000e+04, 1.80000000000000e+04, 0.0, 0.0, 0.0 },
	{ 1.08000000000000e+05, 3.60000000000000e+04, 3.60000000000000e+04, 0.0, 0.0, 0.0 },
	{ 6.34099396490701e+05, 2.11366465496900e+05, 2.11366465496900e+05, 0.0, 0.0, 0.0 },
	{ 1.13477225575052e+06, 3.78257418583506e+05, 3.78257418583506e+05, 0.0, 0.0, 0.0 },
	{ 1.40477225575052e+06, 4.68257418583505e+05, 4.68257418583505e+05, 0.0, 0.0, 0.0 },
	{ 1.67477225575052e+06, 5.58257418583506e+05, 5.58257418583506e+05, 0.0, 0.0, 0.0 },
	{ 1.94477225575052e+06, 6.48257418583506e+05, 6.48257418583506e+05, 0.0, 0.0, 0.0 },
	{ 2.21477225575052e+06, 7.38257418583506e+05, 7.38257418583506e+05, 0.0, 0.0, 0.0 },
	{ 2.48477225575052e+06, 8.28257418583506e+05, 8.28257418583506e+05, 0.0, 0.0, 0.0 }
};


double eps_vs_t(double time, double t_final) {
	const double eps_max = 1.0e-1;
	double eps = eps_max * time;
	return eps;
}


int main(int argc, char **argv)
{

	const int ngp = 1;
	const int n = 2;
	const int dir = 0;
	const double t_final = 0.15;
	const double dt = t_final / time_steps;
	double time = 0.0;

	micropp_params_t mic_params;

	mic_params.ngp = ngp;
	mic_params.size[0] = n;
	mic_params.size[1] = n;
	mic_params.size[2] = n;
	mic_params.type = MIC_HOMOGENEOUS;
	material_set(&mic_params.materials[0], 2, 3.0e7, 0.25, 0.0, 0.0, 1.0e5);
	material_set(&mic_params.materials[1], 0, 3.0e7, 0.25, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[2], 0, 3.0e7, 0.25, 0.0, 0.0, 0.0);
	mic_params.lin_stress = false;

	micropp<3> micro(mic_params);
	micro.print_info();

	ofstream file;
	file.open("result.dat");
	file.precision(14);
	file << scientific;

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

		cout << "sig = ";
		micro.get_stress(0, sig);
		for (int i = 0; i < 6; ++i) {
			assert(fabs(sig_sol[t][i] - sig[i]) < 1.0e-8);
			cout << sig[i] << "\t";
		}
		cout << endl;

		micro.update_vars();

		for (int i = 0; i < 6; ++i) {
			file << eps[i] << "\t";
		}
		for (int i = 0; i < 6; ++i) {
			file << sig[i] << "\t";
		}
		file << endl;

		time += dt;
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "time = " << duration.count() << " ms" << endl;

	return 0;
}
