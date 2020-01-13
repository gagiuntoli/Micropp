/*
 *  This source code is part of MicroPP: a finite element library
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

#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include <cmath>
#include <cassert>
#include <cstring>
#include <ctime>

#ifdef _OPENACC
#include <openacc.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include "util.hpp"
#include "ell.hpp"
#include "material.hpp"
#include "gp.hpp"
#include "instrument.hpp"
#include "params.hpp"
#include "types.hpp"


using namespace std;



template <int tdim>
class micropp {

	protected:
		static constexpr int dim = tdim;                  // 2, 3
		static constexpr int npe = mypow(2, dim);         // 4, 8
		static constexpr int nvoi = dim * (dim + 1) / 2;  // 3, 6
		double bmat_cache[npe][nvoi][npe*dim];

		const int ngp, nx, ny, nz, nn, nndim;
		const int nex, ney, nez, nelem;
		const double lx, ly, lz;
		const double dx, dy, dz;
		const double vol_tot;
		const double wg, ivol, evol;

		const int micro_type, nvars;
		const int nsubiterations;
		const bool subiterations;
		const int mpi_rank;

		gp_t<tdim> *gp_list;

		static const int num_geo_params = 4;
		double geo_params[num_geo_params];

		material_t *material_list[MAX_MATERIALS];
		double ctan_lin_fe[nvoi * nvoi];

		int *elem_type;
		double *elem_stress;
		double *elem_strain;

		const double xg[8][3] = {
			{ -CONSTXG, -CONSTXG, -CONSTXG },
			{ +CONSTXG, -CONSTXG, -CONSTXG },
			{ +CONSTXG, +CONSTXG, -CONSTXG },
			{ -CONSTXG, +CONSTXG, -CONSTXG },
			{ -CONSTXG, -CONSTXG, +CONSTXG },
			{ +CONSTXG, -CONSTXG, +CONSTXG },
			{ +CONSTXG, +CONSTXG, +CONSTXG },
			{ -CONSTXG, +CONSTXG, +CONSTXG } };

		const int nr_max_its;
		const double nr_max_tol;
		const double nr_rel_tol;
		const bool calc_ctan_lin_flag;

		const bool lin_stress;

		/* Linear jacobian for optimization */
		bool use_A0;
		int its_with_A0;
		ell_matrix *A0;

		/* Rule of Mixture Stuff (for 2 mats micro-structure only) */
		double Vm;  // Volume fraction of Matrix
		double Vf;  // Volume fraction of Fiber

		/* IO files */
		const bool write_log_flag;
		int log_id = 0;
		ofstream ofstream_log;

		/* GPU number for device selection */
		int gpu_id = 0;


		/* Private function members */

		/*
		 * Linear homogenizations : Does not perform
		 * FE computation when it is called
		 *
		 * Applies the FE RVE model and Mixture rules
		 * linearly
		 *
		 */
		void homogenize_linear(gp_t<tdim> *gp_ptr);

		/* FE-based homogenizations */
		void homogenize_fe_one_way(gp_t<tdim> *gp_ptr);
		void homogenize_fe_full(gp_t<tdim> *gp_ptr);

		void calc_ctan_lin_fe_models();
		void calc_ctan_lin_mix_rule_Chamis(double ctan[nvoi * nvoi]);

		material_t *get_material(const int e) const;


#pragma acc routine seq
		void get_elem_displ(const double *u,
				    double elem_disp[npe * dim],
				    int ex, int ey, int ez = 0) const;

#pragma acc routine seq
		void get_strain(const double *u, int gp, double strain_gp[nvoi],
				int ex, int ey, int ez = 0) const;

		void get_stress(int gp, const double eps[nvoi],
				const double *vars_old,
				double stress_gp[nvoi],
				int ex, int ey, int ez = 0) const;

		int get_elem_type(int ex, int ey, int ez = 0) const;

		void get_elem_rhs(const double *u, const double *vars_old,
				  double be[npe * dim], int ex, int ey,
				  int ez = 0) const;

		void calc_ave_stress(const double *u, double stress_ave[nvoi],
				     const double *vars_old = nullptr) const;

		void calc_ave_strain(const double *u,
				     double strain_ave[nvoi]) const;

		void calc_fields(double *u, double *vars_old);

		void calc_bmat(int gp, double bmat[nvoi][npe * dim]) const;

		void calc_volume_fractions();

		bool calc_vars_new(const double *u, const double *vars_old,
				   double *vars_new) const;

		newton_t newton_raphson(ell_matrix *A, double *b, double *u,
					double *du, const double strain[nvoi],
					const double *vars_old = nullptr);

		void get_elem_mat(const double *u, const double *vars_old,
				  double Ae[npe * dim * npe * dim],
				  int ex, int ey, int ez = 0) const;

		void set_displ_bc(const double strain[nvoi], double *u);

		double assembly_rhs(const double *u, const double *vars_old,
				    double *b);

		double assembly_rhs_acc(const double *u, const double *vars_old,
					double *b);

		void assembly_mat(ell_matrix *A, const double *u,
				  const double *vars_old);

		// CUDA
		void assembly_mat_cuda(ell_matrix *A, const double *u,
				       const double *vars_old);

		double assembly_rhs_cuda(const double *u, const double *vars_old,
					 double *b);
		//

		void write_vtu(double *u, double *vars_old,
			       const char *filename);

		void write_log();

	public:

		micropp() = delete;

		micropp(const micropp_params_t &params);

		~micropp();

		/* The most important functions */

		void set_strain(const int gp_id, const double *strain);

		void get_stress(const int gp_id, double *stress) const;

		void get_ctan(const int gp_id, double *ctan) const;

		void homogenize();

		void homogenize_linear();

		/* Extras */

		int is_non_linear(const int gp_id) const;

		int get_non_linear_gps(void) const;

		int get_cost(int gp_id) const;

		bool has_converged(int gp_id) const;

		bool has_subiterated(int gp_id) const;

		void output(int gp_id, const char *filename);

		void output2(const int gp_id, const int elem_global,
			     const int time_step);

		void update_vars();

		void write_restart(const int restart_id) const;

		void read_restart(const int restart_id) const;

		void print_info() const;

};
