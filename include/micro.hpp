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

#ifndef MICRO_HPP
#define MICRO_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <map>

#include <cmath>
#include <cassert>
#include <cstring>
#include <ctime>

#include "util.hpp"
#include "ell.hpp"
#include "material.hpp"
#include "gp.hpp"
#include "instrument.hpp"

#ifdef _OPENACC
#include <openacc.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_DIM         3
#define NUM_VAR_GP      7  // eps_p_1 (6) , alpha_1 (1)
#define MAX_MATERIALS   3

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

	void print()
	{
		cout << "newton.its        : " << its << endl;
		cout << "newton.solver_its : " << solver_its << endl;
		cout << "newton.converged  : " << converged << endl;
	}

} newton_t;


typedef struct {

	int ngp = 1;
	int size[3];
	int type = 0;
	double geo_params[4] = {0.1, 0.1, 0.1, 0.1};
	struct material_base materials[4];
	int *coupling = nullptr;
	bool subiterations = false;
	int nsubiterations = 10;
	int mpi_rank = 0;
	int nr_max_its = NR_MAX_ITS;
	double nr_max_tol = NR_MAX_TOL;
	double nr_rel_tol = NR_REL_TOL;
	int cg_max_its = CG_MAX_ITS;
	double cg_abs_tol = CG_ABS_TOL;
	double cg_rel_tol = CG_REL_TOL;
	bool calc_ctan_lin = true;
	bool use_A0 = false;
	int its_with_A0 = 1;
	bool lin_stress = true;
	bool write_log = false;

	void print()
	{
		cout << "ngp  : " << ngp << endl;
		cout << "size : " << size[0] << endl;
		cout << "type  : " << type << endl;
		cout << "geo_params : " << geo_params[0] << endl;
		cout << "subiterations : " << subiterations << endl;
		cout << "nsubiterations : " << nsubiterations << endl;
		cout << "mpi_rank : " << mpi_rank << endl;
		cout << "nr_max_its : " << nr_max_its << endl;
		cout << "nr_max_tol : " << nr_max_tol << endl;
		cout << "nr_rel_tol : " << nr_rel_tol << endl;
		cout << "calc_ctan_lin : " << calc_ctan_lin << endl;
		cout << "use_A0 : " << use_A0 << endl;
		cout << "its_with_A0 : " << its_with_A0 << endl;
		cout << "lin_stress : " << lin_stress << endl;
		cout << "write_log : " << write_log << endl;
	}

} micropp_params_t;


enum {
	MIC_HOMOGENEOUS,
	MIC_SPHERE,
	MIC_LAYER_Y,
	MIC_CILI_FIB_X,
	MIC_CILI_FIB_Z,
	MIC_CILI_FIB_XZ,
	MIC_QUAD_FIB_XYZ,
	MIC_QUAD_FIB_XZ,
	MIC_QUAD_FIB_XZ_BROKEN_X,
	MIC3D_SPHERES,
	MIC3D_8,
	MIC3D_FIBS_20_DISORDER
};

static map<int, std::string> micro_names = {
	{MIC_HOMOGENEOUS, "MIC_HOMOGENEOUS"},
	{MIC_SPHERE, "MIC_SPHERE"},
	{MIC_LAYER_Y, "MIC_LAYER_Y"},
	{MIC_CILI_FIB_X, "MIC_CILI_FIB_X"},
	{MIC_CILI_FIB_Z, "MIC_CILI_FIB_Z"},
	{MIC_CILI_FIB_XZ, "MIC_CILI_FIB_XZ"},
	{MIC_QUAD_FIB_XYZ, "MIC_QUAD_FIB_XYZ"},
	{MIC_QUAD_FIB_XZ, "MIC_QUAD_FIB_XZ"},
	{MIC_QUAD_FIB_XZ_BROKEN_X, "MIC_QUAD_FIB_XZ_BROKEN_X"},
	{MIC3D_SPHERES, "MIC3D_SPHERES"},
	{MIC3D_8, "MIC3D_8"},
	{MIC3D_FIBS_20_DISORDER, "MIC3D_FIBS_20_DISORDER"}
};

/*
 * MIC_SPHERES : (2 materials) One sphere in the middle
 *
 * MIC_HOMOGENEOUS : Only one material (mat[0])
 *
 * MIC3D_SPHERES : (2 materials) Random spheres.
 *
 * MIC3D_8 : (3 materiales) 2 cilinders at 90 deg with a layer around the
 * perimeter and a flat layer between the fibers.
 *
 * MIC3D_FIBS_20_DISORDER: (2 materiales) 20 fibers in random directions.
 *
 */


enum {
	FE_LINEAR,
	FE_ONE_WAY,
	FE_FULL,
	MIX_RULE_CHAMIS
};


static map<int, int> gp_counter = {
	{FE_LINEAR, 0},
	{FE_ONE_WAY, 0},
	{FE_FULL, 0},
	{MIX_RULE_CHAMIS, 0}
};


using namespace std;


template <int tdim>
class micropp {

	protected:
		static constexpr int dim = tdim;                  // 2, 3
		static constexpr int npe = mypow(2, dim);         // 4, 8
		static constexpr int nvoi = dim * (dim + 1) / 2;  // 3, 6
		double calc_bmat_cache[npe][nvoi][npe*dim];

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
		 * Linear homogenizations 
		 * Applies to FE RVE model and Mixture rules
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
		void get_elem_nodes(int n[npe],
				    int ex, int ey, int ez = 0) const;

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

		void assembly_mat_acc(ell_matrix *A, const double *u,
				      const double *vars_old);

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


#endif
