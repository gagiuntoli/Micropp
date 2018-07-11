/*
 * This source code is part of MicroPP: a finite element library
 * to solve microstructural problems for composite materials.
 *
 * Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                        Guido Giuntoli <gagiuntoli@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>

#include <cmath>

#include "ell.hpp"

#define MAX_DIM       3
#define MAX_MAT_PARAM 10
#define MAX_MATS      10
#define MAX_GP_VARS   10
#define INT_VARS_GP   7		// eps_p_1, alpha_1
#define NUM_VAR_GP    7		// eps_p_1, alpha_1

#define CG_MAX_TOL    1.0e-8
#define CG_MAX_ITS    2000
#define NR_MAX_TOL    1.0e-5
#define NR_MAX_ITS    40

#define glo_elem3D(ex,ey,ez) ((ez) * (nx-1) * (ny-1) + (ey) * (nx-1) + (ex))
#define intvar_ix(e,gp,var) ((e) * 8 * INT_VARS_GP + (gp) * INT_VARS_GP + (var))

using namespace std;

struct gp_t {
	int id;
	int nr_its[7];
	double *int_vars_n;
	double *int_vars_k;
	double macro_strain[6];
	double macro_stress[6];
	double macro_ctan[36];
	double nr_err[7];
	double inv_max;
};

struct material_t {
	double E;
	double nu;
	double k;
	double mu;
	double lambda;
	double Ka;
	double Sy;
	bool plasticity;
	bool damage;
};

class micropp_t {

	private:
		const int dim;
		const int nx, ny, nz, nn;
		const int nex, ney, nez;
		const double lx, ly, lz, dx, dy, dz, width, inv_tol;
		const int npe, nvoi, nelem;
		const int micro_type, num_int_vars;

		bool output_files_header;

		double micro_params[5];
		int numMaterials;
		material_t material_list[MAX_MATS];
		double ctan_lin[36];

		list<gp_t> gauss_list;

		ell_matrix A;
		double * u;
		double * du;
		double * b;
		ell_solver solver;

		double * elem_stress;
		double * elem_strain;
		int * elem_type;
		double * vars_old;
		double * vars_new;

		double inv_max;

	public:
		micropp_t(const int dim, const int size[3], const int micro_type,
		          const double *micro_params, const int *mat_types,
		          const double *params);

		~micropp_t();

		void calc_ctan_lin();

		bool is_linear(const double *macro_strain);

		double get_inv_1(const double *tensor);
		double get_inv_2(const double *tensor);

		void set_macro_strain(const int gp_id, const double *macro_strain);
		void get_macro_stress(const int gp_id, double *macro_stress);
		void get_macro_ctan(const int gp_id, double *macro_ctan);

		void homogenize();
		void update_vars();
		void get_nl_flag(int gp_id, int *nl_flag);

		void set_displ(double *eps);
		double assembly_rhs(bool *nl_flag);

		void get_elem_rhs2D(int ex, int ey, bool *nl_flag, double (&be)[2 * 4]);
		void get_elem_rhs3D(int ex, int ey, int ez, bool *nl_flag, double (&be)[3 * 8]);

		void assembly_mat();

		void get_elem_mat2D(int ex, int ey, double (&Ae)[2 * 4 * 2 * 4]);
		void get_elem_mat3D(int ex, int ey, int ez, double (&Ae)[3 * 8 * 3 * 8]);

		void solve();
		void newton_raphson(bool *nl_flag, int *its, double *err);

		void get_ctan_plast_sec(int ex, int ey, int ez, int gp, double ctan[6][6]);
		void get_ctan_plast_exact(int ex, int ey, int ez, int gp, double ctan[6][6]);
		void get_ctan_plast_pert(int ex, int ey, int ez, int gp, double ctan[6][6]);

		void get_strain2D(int ex, int ey, int gp, double *strain_gp);
		void get_strain3D(int ex, int ey, int ez, int gp, double *strain_gp);

		void get_stress2D(int ex, int ey, int gp, double strain_gp[3],
		                  bool *nl_flag, double *stress_gp);

		void get_stress3D(int ex, int ey, int ez, int gp, double strain_gp[3],
		                  bool *nl_flag, double *stress_gp);

		void get_dev_tensor(double tensor[6], double tensor_dev[6]);
		void plastic_step(const material_t *material, double eps[6],
		                  double eps_p_1[6], double alpha_1,
		                  double eps_p[6], double *alpha,
		                  bool *nl_flag, double stress[6]);

		void getElemDisp(int ex, int ey, double *elem_disp);
		void getElemDisp(int ex, int ey, int ez, double *elem_disp);

		int get_elem_type2D(int ex, int ey);
		int get_elem_type3D(int ex, int ey, int ez);

		material_t get_material(const int e);

		void calc_bmat_3D(int gp, double bmat[6][3 *8]);

		void calc_fields();
		void calc_ave_stress(double stress_ave[6]);
		void calc_ave_strain(double strain_ave[6]);

		void output(int tstep, int gp_id);
		void write_vtu(int tstep, int gp_id);
		void write_info_files();
};
