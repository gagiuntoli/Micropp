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

#ifndef MICRO_HPP
#define MICRO_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>
#include <cassert>
#include <cstring>
#include <ctime>

#include "util.hpp"
#include "ell.hpp"
#include "material.hpp"
#include "gp.hpp"
#include "instrument.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_DIM         3
#define MAX_MATS        10
#define NUM_VAR_GP      7  // eps_p_1 (6) , alpha_1 (1)

#define CG_MIN_ERR      1.0e-50
#define CG_MAX_ITS      1000
#define CG_REL_ERR      1.0e-5

#define FILTER_REL_TOL  1.0e-5

#define D_EPS_CTAN_AVE  1.0e-8

#define CONSTXG         0.577350269189626

#define glo_elem(ex,ey,ez)   ((ez) * (nx-1) * (ny-1) + (ey) * (nx-1) + (ex))
#define intvar_ix(e,gp,var)  ((e) * npe * NUM_VAR_GP + (gp) * NUM_VAR_GP + (var))

#define NR_MAX_TOL      1.0e-10
#define NR_MAX_ITS      4
#define NR_REL_TOL      1.0e-3 // factor against first residual


typedef struct {

	/* Results from newton-raphson loop */
	int its = 0;
	int solver_its = 0;
	bool converged = false;

} newton_t;


enum {
       	MIC_SPHERE,
       	MIC_LAYER_Y,
       	MIC_CILI_FIB_Z,
       	MIC_CILI_FIB_XZ,
       	MIC_QUAD_FIB_XYZ,
       	MIC_QUAD_FIB_XZ,
	MIC_QUAD_FIB_XZ_BROKEN_X,
	MIC_SPHERES
};

enum {
       	NO_COUPLING,
       	ONE_WAY,
       	FULL
};


using namespace std;


template <int tdim>
class micropp {

	protected:
		static constexpr int dim = tdim;                  // 2, 3
		static constexpr int npe = mypow(2, dim);         // 4, 8
		static constexpr int nvoi = dim * (dim + 1) / 2;  // 3, 6

		const int ngp, nx, ny, nz, nn, nndim;
		const int nex, ney, nez, nelem;
		const double lx, ly, lz;
		const double dx, dy, dz;
		const double vol_tot;
		const double special_param, wg, ivol;

		const int micro_type, nvars;
		const int nsubiterations;
		const bool subiterations;

		gp_t<tdim> *gp_list;

		int coupling;

		double micro_params[5];
		int numMaterials;
		material_t *material_list[MAX_MATS];
		double ctan_lin[nvoi * nvoi];

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

		void calc_ctan_lin();
		material_t *get_material(const int e) const;

		void get_elem_nodes(int n[npe],
				    int ex, int ey, int ez = 0) const;

		void get_elem_displ(const double *u,
				    double elem_disp[npe * dim],
				    int ex, int ey, int ez = 0) const;

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

		void assembly_mat(ell_matrix *A, const double *u,
				  const double *vars_old);

		void write_vtu(double *u, double *vars_old,
			       const char *filename);

	public:

		micropp() = delete;

		micropp(const int ngp, const int size[3], const int micro_type,
			const double *micro_params, const material_t *materials,
			const int _coupling = ONE_WAY,
			const bool _subiterations = false,
			const int _nsubiterations = 10,
			const int max_its = NR_MAX_ITS,
			const double max_tol = NR_MAX_TOL,
			const double rel_tol = NR_REL_TOL);

		~micropp();

		/* The most important functions */

		void set_strain(const int gp_id, const double *strain);

		void get_stress(const int gp_id, double *stress) const;

		void get_ctan(const int gp_id, double *ctan) const;

		void homogenize();

		/* Extras */

		int is_non_linear(const int gp_id) const;

		int get_non_linear_gps(void) const;

		int get_cost(int gp_id) const;

		bool has_converged(int gp_id) const;

		bool has_subiterated(int gp_id) const;

		void output(int gp_id, const char *filename);

		void update_vars();

		void print_info() const;

};


#endif
