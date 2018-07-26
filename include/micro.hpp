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

#include <cmath>
#include <cassert>
#include <cstring>

#include "util.hpp"
#include "ell.hpp"
#include "gp.hpp"
#include "instrument.hpp"

#define MAX_DIM         3
#define MAX_MAT_PARAM   10
#define MAX_MATS        10
#define MAX_GP_VARS     10
#define INT_VARS_GP     7  // [eps_p_1(gp1), alpha_1(gp1), eps_p_1(gp2), ...]
#define NUM_VAR_GP      7  // eps_p_1 (6) , alpha_1 (1)

#define CG_MIN_ERR      1.0e-8
#define CG_MAX_ITS      2000
#define NR_MAX_TOL      1.0e-5
#define NR_MAX_ITS      40
#define D_EPS_CTAN      1.0e-8
#define D_EPS_CTAN_AVE  1.0e-8

#define CONSTXG         0.577350269189626

#define glo_elem(ex,ey,ez)   ((ez) * (nx-1) * (ny-1) + (ey) * (nx-1) + (ex))
#define intvar_ix(e,gp,var)  ((e) * 8 * INT_VARS_GP + (gp) * INT_VARS_GP + (var))

using namespace std;

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


template <int tdim>
class micropp {

	private:
		// Variables (static constexpr)
		static constexpr int dim = tdim;
		static constexpr int npe = mypow(2, dim);         // 4, 8
		static constexpr int nvoi = dim * (dim + 1) / 2;  // 3, 6

		// Constants only vars
		const int ngp, nx, ny, nz, nn;
		const int nex, ney, nez, nelem;
		const double lx, ly, lz;
		const double dx, dy, dz, dxi, dyi, dzi;
		const double width, inv_tol;

		const int micro_type, num_int_vars;
		gp_t *gp_list;

		// Other variables
		bool output_files_header;

		double micro_params[5];
		int numMaterials;
		material_t material_list[MAX_MATS];
		double ctan_lin[36];

		ell_matrix A;

		double *b;
		double *du;
		double *u;
		double *u_aux;

		int *elem_type;
		double *elem_stress;
		double *elem_strain;
		double *vars_old_aux, *vars_old;
		double *vars_new_aux, *vars_new;

		double inv_max;

		const double xg[8][3] = { {-CONSTXG, -CONSTXG, -CONSTXG},
		                          {+CONSTXG, -CONSTXG, -CONSTXG},
		                          {+CONSTXG, +CONSTXG, -CONSTXG},
		                          {-CONSTXG, +CONSTXG, -CONSTXG},
		                          {-CONSTXG, -CONSTXG, +CONSTXG},
		                          {+CONSTXG, -CONSTXG, +CONSTXG},
		                          {+CONSTXG, +CONSTXG, +CONSTXG},
		                          {-CONSTXG, +CONSTXG, +CONSTXG} };

	protected:
		// Common
		void calc_ctan_lin();
		bool is_linear(const double *macro_strain);
		double get_inv_1(const double *tensor) const;
		material_t get_material(const int e) const;

		void get_strain(int gp, double *strain_gp,
				int ex, int ey, int ez = 0) const;

		void get_elem_nodes(int n[8], int ex, int ey, int ez = 0) const;

		void get_elem_displ(const double *u, double *elem_disp,
				int ex, int ey, int ez = 0) const;

		void isolin_get_stress(
				const material_t *material, const double eps[6],
				double stress[6]) const;

		void plastic_get_stress(
				const material_t *material, const double eps[6],
				const double eps_p_old[6], double alpha_old,
				double stress[6]) const;

		void get_stress(int gp, const double eps[nvoi],
				double stress_gp[nvoi], int ex, int ey, int ez = 0) const;

		// Specialized
		template <typename... Rest>
		int get_elem_type(Rest...);

		template <typename... Rest>
		void get_elem_rhs(double *be, Rest...) const;

		template <typename... Rest>
		void get_elem_mat(Rest...) const;

		void set_displ_bc(const double *eps);
		double assembly_rhs();
		void assembly_mat();

		template <typename T>
		void calc_bmat(int gp, T bmat) const;

		void newton_raphson(int *_its, double *_err);

		void calc_ave_stress(double stress_ave[nvoi]) const;
		void calc_ave_strain(double strain_ave[nvoi]) const;

		bool calc_vars_new();

		void calc_fields();

		void write_vtu(int tstep, int gp_id);

		// Functions Only for 3D

		void get_dev_tensor(const double tensor[6], double tensor_dev[6]) const;

		bool plastic_law(
				const material_t *material, const double eps[6],
				const double eps_p_old[6], double alpha_old,
				double *_dl, double _normal[6], double _s_trial[6]) const;

		void plastic_get_ctan(
				const material_t *material, const double eps[6],
				const double eps_p_old[6], double alpha_old,
				double ctan[6][6]) const;

		bool plastic_evolute(
				const material_t *material, const double eps[6],
				const double eps_p_old[6], double alpha_old, 
				double *eps_p_new, double *alpha_new) const;

		void isolin_get_ctan(
				const material_t *material, double ctan[6][6]) const;

		void initialize(const double *micro_params, const int *mat_types,
				const double *params);

	public:
		micropp(const int ngp,const int size[3], const int micro_type,
		        const double *micro_params, const int *mat_types,
		        const double *params);
		~micropp();

		// common Functions

		int get_nl_flag(const int gp_id) const;
		void set_macro_strain(const int gp_id, const double *macro_strain);
		void get_macro_stress(const int gp_id, double *macro_stress) const;
		void get_macro_ctan(const int gp_id, double *macro_ctan) const;
		void homogenize();
		void output(int tstep, int gp_id);
		void write_info_files();
		void update_vars();
		void print_info() const;
};
