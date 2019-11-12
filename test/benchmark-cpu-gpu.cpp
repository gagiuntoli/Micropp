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
#include <chrono>


#include "micropp.hpp"


using namespace std;
using namespace std::chrono;


const double strain[6] = { 0.001, 0.0, 0., 0., 0., 0. };

class test_t : public micropp<3> {

	public:
		test_t(micropp_params_t mic_params):micropp<3>(mic_params) {};

		~test_t() {};

		void newton_raphson(ofstream & file)
		{
			double cg_err;

			ell_matrix A;  // Jacobian
			const int ns[3] = { nx, ny, nz };
			ell_init(&A, dim, dim, ns, CG_ABS_TOL, CG_REL_TOL, CG_MAX_ITS);
			double *b = (double *) calloc(nndim, sizeof(double));
			double *du = (double *) calloc(nndim, sizeof(double));
			double *u = (double *) calloc(nndim, sizeof(double));

			cout <<  "CPU Execution" << endl;
			memset(u, 0.0, nndim * sizeof(double));

			set_displ_bc(strain, u);

			auto time_1 = high_resolution_clock::now();

			double norm = assembly_rhs_acc(u, nullptr, b);

			auto time_2 = high_resolution_clock::now();

			cout << "|r| : " << norm << endl;

			auto time_3 = high_resolution_clock::now();

			assembly_mat(&A, u, nullptr);

			auto time_4 = high_resolution_clock::now();
			auto time_5 = high_resolution_clock::now();

			int cg_its = ell_solve_cgpd(&A, b, du, &cg_err);

			auto time_6 = high_resolution_clock::now();

			cout << "CG Its : " << cg_its << endl;
			cout << "CG Err : " << cg_err << endl;

			for (int i = 0; i < nn * dim; ++i)
				u[i] += du[i];

			auto time_7 = high_resolution_clock::now();

			norm = assembly_rhs_acc(u, nullptr, b);

			auto time_8 = high_resolution_clock::now();

			cout << "|r| : " << norm << endl;

			auto ass_res = 
				duration_cast<milliseconds>(time_2 - time_1) +
				duration_cast<milliseconds>(time_8 - time_7);

			auto ass_mat = 
				duration_cast<milliseconds>(time_4 - time_3);

			auto solver = 
				duration_cast<milliseconds>(time_6 - time_5);

			auto ass_tot = ass_res.count() + ass_mat.count();
			auto total = solver.count() + ass_tot;
			double percentage_ass = (100.0 * ass_tot) / total;
			double percentage_sol = (100.0 * solver.count()) / total;

			cout << "ass_res : " << ass_res.count() << " ms" << endl;
			cout << "ass_mat : " << ass_mat.count() << " ms" << endl;
			cout << "ass_tot : " << ass_tot << " ms" << endl;
			cout << "solver : " << solver.count() << " ms" << endl;
			cout
				<< "ass : " << percentage_ass << " \% " 
				<< "sol : " << percentage_sol << " \%" << endl;

			file
				<< nx - 1 << "\t"
				<< ass_res.count() + ass_mat.count() << "\t"
				<< solver.count() << "\t"
				<< total << "\t"
				<< percentage_ass << "\t"
				<< percentage_sol << "\t";

		//------------------------------------------------------------	
		
			cout <<  "GPU Execution" << endl;
			memset(u, 0.0, nndim * sizeof(double));

			set_displ_bc(strain, u);

			time_1 = high_resolution_clock::now();

			norm = assembly_rhs_acc(u, nullptr, b);

			time_2 = high_resolution_clock::now();

			cout << "|r| : " << norm << endl;

			time_3 = high_resolution_clock::now();

			assembly_mat_acc(&A, u, nullptr);

			time_4 = high_resolution_clock::now();
			time_5 = high_resolution_clock::now();

			cg_its = ell_solve_cgpd_acc(&A, b, du, &cg_err);

			time_6 = high_resolution_clock::now();

			cout << "CG Its : " << cg_its << endl;
			cout << "CG Err : " << cg_err << endl;

			for (int i = 0; i < nn * dim; ++i)
				u[i] += du[i];

			time_7 = high_resolution_clock::now();

			norm = assembly_rhs_acc(u, nullptr, b);

			time_8 = high_resolution_clock::now();

			cout << "|r| : " << norm << endl;

			ass_res = 
				duration_cast<milliseconds>(time_2 - time_1) +
				duration_cast<milliseconds>(time_8 - time_7);

			ass_mat = 
				duration_cast<milliseconds>(time_4 - time_3);

			solver = 
				duration_cast<milliseconds>(time_6 - time_5);

			ass_tot = ass_res.count() + ass_mat.count();
			total = solver.count() + ass_tot;
			percentage_ass = (100.0 * ass_tot) / total;
			percentage_sol = (100.0 * solver.count()) / total;

			cout << "ass_res : " << ass_res.count() << " ms" << endl;
			cout << "ass_mat : " << ass_mat.count() << " ms" << endl;
			cout << "ass_tot : " << ass_tot << " ms" << endl;
			cout << "solver : " << solver.count() << " ms" << endl;
			cout
				<< "ass : " << percentage_ass << " \% " 
				<< "sol : " << percentage_sol << " \%" << endl;

			file
				<< ass_res.count() + ass_mat.count() << "\t"
				<< solver.count() << "\t"
				<< total << "\t"
				<< percentage_ass << "\t"
				<< percentage_sol << endl;

			ell_free(&A);
			free(b);
			free(u);
			free(du);
		};

};


int main (int argc, char *argv[])
{

#ifndef _OPENACC
	cout << "Warning : The GPU code will be executed all in CPU because "
		"the code wasn't compiled with OpenACC." << endl;
#endif

	micropp_params_t mic_params;
	mic_params.ngp = 1;
	mic_params.type = MIC3D_8;
	mic_params.geo_params[0] = 0.1;
	mic_params.geo_params[1] = 0.02;
	mic_params.geo_params[2] = 0.01;
	material_set(&mic_params.materials[0], 0, 3.0e7, 0.25, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[1], 0, 3.0e8, 0.25, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[2], 0, 3.0e7, 0.25, 0.0, 0.0, 0.0);
	mic_params.mpi_rank = 0;
	mic_params.calc_ctan_lin = false;
	mic_params.lin_stress = false;

	const int n_cases = 10;
	const int n[n_cases] = { 11, 21, 31, 41, 51, 61, 71, 81, 91, 101 };

	ofstream file;
	file.open("benchmark-cpu-gpu.dat");
	file << "#             CPU                             GPU" << endl;
	file 
		<< "#N   time_ass  time_sol  total   ass\%   sol\% "
		<< "time_ass  time_sol  total   ass\%   sol\%" << endl;

	for (int i = 0; i < n_cases; ++i) {
		mic_params.size[0] = n[i];
		mic_params.size[1] = n[i];
		mic_params.size[2] = n[i];

		test_t test(mic_params);
		test.newton_raphson(file);
	}

	return 0;
}
