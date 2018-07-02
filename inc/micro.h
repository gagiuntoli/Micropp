/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
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

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>
#include <cmath>
#include "ell.h"

#define MAX_MAT_PARAM 10
#define MAX_MATS      10
#define MAX_GP_VARS   10
#define INT_VARS_GP   7  // eps_p_1, alpha_1
#define NUM_VAR_GP    7  // eps_p_1, alpha_1

#define glo_elem3D(ex,ey,ez) ((ez)*(nx-1)*(ny-1) + (ey)*(nx-1) + (ex))
#define intvar_ix(e,gp,var) ((e)*8*INT_VARS_GP + (gp)*INT_VARS_GP + (var))

using namespace std;

struct convergence_t {
  	int NR_Its_Stress;
  	double NR_Err_Stress;
  	int NR_Its_Ctan[6];
  	double NR_Err_Ctan[6];
  	double I_reached, I_reached_aux;
};

struct GaussPoint_t {
  	int id;
  	double *int_vars_n;
  	double *int_vars_k;
  	double macro_strain[6];
  	double macro_stress[6];
  	double macro_ctan[36];
  	convergence_t convergence;
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

class Problem {

  	private:

    	int dim, npe, nvoi;
    	int nx, ny, nz, nn; 
    	double lx, ly, lz, dx, dy, dz;
    	int nelem;
    	int size_tot;

    	int micro_type;
    	double micro_params[5];
    	int numMaterials;
    	material_t material_list[MAX_MATS];
    	std::list<GaussPoint_t> gauss_list;
    	double ctan_lin[6][6];

    	double *elem_strain;
    	double *elem_stress;
    	int *elem_type;
    	int num_int_vars;
    	double *vars_old, *vars_new; 

    	ell_matrix A;
    	double *u, *du, *b;

    	double I_max;
    	double I_reached;

    	double NR_norm;
    	int NR_its, NR_non_linear;

  	public:

    	Problem (int dim, int size[3], int micro_type, double *micro_params, int *mat_types, double *params);
    	~Problem();

    	void calcCtanLinear();

    	bool LinearCriteria (const double *macro_strain);
    	double Invariant_I1 (const double *tensor);
    	double Invariant_I2 (const double *tensor);

    	void set_macro_strain(const int gp_id, const double *macro_strain);
    	void get_macro_stress(const int gp_id, double *macro_stress);
    	void get_macro_ctan(const int gp_id, double *macro_ctan);
    	void homogenize();
    	void update_vars();

    	void solve();
    	void newtonRaphson (bool *non_linear_flag);

    	void getNonLinearFlag (int macroGp_id, int *non_linear);
    	void getIntVars (int macroGp_id, int n, int *int_vars);

    	void setDisp (double *eps);
    	void Assembly_A ();
    	double Assembly_b (bool *non_linear_flag);

    	void getElemental_A (int ex, int ey, double (&Ae)[2*4*2*4]);
    	void getElemental_A (int ex, int ey, int ez, double (&Ae)[3*8*3*8]);
    	void getCtanPlasSecant (int ex, int ey, int ez, int gp, double ctan[6][6]);
    	void getCtanPlasExact (int ex, int ey, int ez, int gp, double ctan[6][6]);
    	void getCtanPlasPert (int ex, int ey, int ez, int gp, double ctan[6][6]);

    	void getElemental_b (int ex, int ey, bool *non_linear_flag, double (&be)[2*4]);
    	void getElemental_b (int ex, int ey, int ez, bool *non_linear_flag, double (&be)[3*8]);

    	void getStrain (int ex, int ey, int gp, double *strain_gp);
    	void getStrain (int ex, int ey, int ez, int gp, double *strain_gp);

    	void getStress (int ex, int ey, int gp, double strain_gp[3], bool *non_linear_flag, double *stress_gp);
    	void getStress (int ex, int ey, int ez, int gp, double strain_gp[3], bool *non_linear_flag, double *stress_gp);
    	void getDeviatoric (double tensor[6], double tensor_dev[6]);
    	void plastic_step(
				material_t &material, double eps[6], double eps_p_1[6], double alpha_1, double eps_p[6], 
				double *alpha, bool *non_linear, double stress[6]);

    	void getElemDisp (int ex, int ey, double *elem_disp);
    	void getElemDisp (int ex, int ey, int ez, double *elem_disp);

    	int getElemType (int ex, int ey);
    	int getElemType (int ex, int ey, int ez);
    	void getMaterial (int e, material_t &material);

    	void calc_bmat_3D (int gp, double bmat[6][3*8]);

    	void calcDistributions ();
    	void calcAverageStress (double stress_ave[6]);
    	void calcAverageStrain (double strain_ave[6]);

    	void output (int time_step, int Gauss_ID, double *macro_strain);
    	void writeVtu (int time_step, int elem);
    	void writeConvergenceFile ();
    	bool output_files_header;
};
