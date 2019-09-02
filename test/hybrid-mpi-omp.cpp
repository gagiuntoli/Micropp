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
#include <mpi.h>


#include "micropp.hpp"


using namespace std;
using namespace std::chrono;


#define D_EPS 5.0e-4


int main(int argc, char **argv)
{

	int rank, nproc;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " n [ngp] [steps]" << endl;
		return(1);
	}

	const double perc_nol = 0.25;
	const int dir = 1;
	const int n = atoi(argv[1]);
	const int ngp = (argc > 2 ? atoi(argv[2]) : 10);
	const int time_steps = (argc > 3 ? atoi(argv[3]) : 10);  // Optional value
	const int ngp_nol = perc_nol * ngp;
	const int ngp_lin = ngp - ngp_nol;
	if (!rank) {
		cout << "ngp_nol = " << ngp_nol << " ngp_lin = " << ngp_lin << endl;
	}

	assert(n > 1 && ngp > 0 && time_steps >= 0);
	const int ngp_per_mpi = ngp / nproc + ((ngp % nproc > rank) ? 1 : 0);
	const int gp_start = (ngp % nproc > rank) ? (ngp / nproc + 1) * rank : (ngp / nproc + 1) * (ngp % nproc) + (ngp / nproc) * (rank - ngp % nproc);
	const int gp_end = gp_start + ngp_per_mpi;

	const int ngp_nol_loc = (gp_start > ngp_nol) ? 0 : ((gp_end >= ngp_nol) ? ngp_nol - gp_start : ngp_nol);
	const int ngp_lin_loc = ngp_per_mpi - ngp_nol_loc;

	cout    << "RANK = " << rank << " ngp = " << ngp_per_mpi 
	        << " gp_start = " << gp_start << " gp_end = " << gp_end
		<< " ngp_nol_loc = " << ngp_nol_loc << " ngp_lin_loc = " << ngp_lin_loc << endl;

	micropp_params_t mic_params;

	mic_params.ngp = ngp_per_mpi;
	mic_params.size[0] = n;
	mic_params.size[1] = n;
	mic_params.size[2] = n;
	mic_params.type = MIC_SPHERE;
	material_set(&mic_params.materials[0], 0, 1.0e6, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[1], 2, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[2], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	mic_params.mpi_rank = rank;
	mic_params.calc_ctan_lin = false;
	mic_params.lin_stress = false;

	//mic_params.print();

	micropp<3> micro(mic_params);
	//micro.print_info();

	auto start = high_resolution_clock::now();

	double sig[6];
	double eps_lin[6] = { 0. };
	double eps_nol[6] = { 0. };

	cout << scientific;

	for (int t = 0; t < time_steps; ++t) {

		cout << "Time step = " << t << endl;

		eps_nol[dir] += D_EPS * 10;
		eps_lin[dir] += D_EPS * 0.1;

		for (int gp = 0; gp < ngp_per_mpi; ++gp) {
			if (gp < ngp_nol_loc) {
				micro.set_strain(gp, eps_nol);
			} else {
				micro.set_strain(gp, eps_lin);
			}
		}
		cout << "eps_lin = ";
		for (int i = 0; i < 6; ++i) {
			cout << eps_lin[i] << " ";
		}
		cout << endl;
		cout << "eps_nol = ";
		for (int i = 0; i < 6; ++i) {
			cout << eps_nol[i] << " ";
		}
		cout << endl;

		cout << "Homogenizing ..." << endl;
		micro.homogenize();

		for (int gp = 0; gp < ngp_per_mpi; ++gp) {
			micro.get_stress(gp, sig);
			int non_linear = micro.is_non_linear(gp);
			int cost = micro.get_cost(gp);
			bool has_converged = micro.has_converged(gp);
			cout << "NL = " << non_linear << " ";
			cout << "sig = ";
			for (int i = 0; i < 6; ++i) {
				cout << sig[i] << "\t";
			}
			cout << "cost = " << cost << " ";
			cout << "conv = " << has_converged << " ";
			cout << endl;
		}
		cout << endl;

		micro.update_vars();

		MPI_Barrier(MPI_COMM_WORLD);
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "time = " << duration.count() << " ms" << endl;

	MPI_Finalize();

	return 0;
}
