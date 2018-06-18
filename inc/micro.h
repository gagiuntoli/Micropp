#include <vector>
#include <iostream>
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

struct MacroGp_t {
  int id;
  bool non_linear;
  bool non_linear_aux;
  double *int_vars;
  double *int_vars_aux;
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

  public:

    int dim, npe, nvoi;
    int nx, ny, nz, nn; 
    double lx, ly, lz, dx, dy, dz;
    int nelem;
    int size_tot;

    int micro_type;
    double micro_params[5];
    int numMaterials;
    material_t material_list[MAX_MATS];
    std::list<MacroGp_t> MacroGp_list;
    double CtanLinear[6][6];

    double *elem_strain;
    double *elem_stress;
    int *elem_type;
    int num_int_vars;
    double *vars_old, *vars_new; 

    ell_matrix A;
    double *u, *du, *b;

    Problem (int dim, int size[3], int micro_type, double *micro_params, int *mat_types, double *params);
    ~Problem (void);
    void calcCtanLinear (void);

    int LinCriteria;
    void loc_hom_Stress (int macro_id, double *MacroStrain, double *MacroStress);
    void loc_hom_Ctan (int macroGp_id, double *MacroStrain, double *MacroCtan);
    bool LinearCriteria (double *MacroStrain);
    double Invariant_I1 (double *tensor);
    double Invariant_I2 (double *tensor);
    void updateIntVars (void);

    double NR_norm;
    int NR_its, NR_non_linear;
    void solve (void);
    void newtonRaphson (bool *non_linear_flag);

    void getParams_LinCriteria (int *LinCriteria) {*LinCriteria = this->LinCriteria;};
    void getParams_NR (int *NR_its, double *NR_norm, int *NR_non_linear ) 
    {*NR_its = this->NR_its; *NR_norm = this->NR_norm; *NR_non_linear = this->NR_non_linear;};
    void getNonLinearFlag (int macroGp_id, int *non_linear);
    void getIntVars (int macroGp_id, int n, int *int_vars);

    void setDisp (double *eps);

    void Assembly_A (void);
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
    void plasticStep(
	material_t &material, double eps[6], double eps_p_1[6], double alpha_1, double eps_p[6], 
	double *alpha, bool *non_linear, double stress[6]);

    void getElemDisp (int ex, int ey, double *elem_disp);
    void getElemDisp (int ex, int ey, int ez, double *elem_disp);

    int getElemType (int ex, int ey);
    int getElemType (int ex, int ey, int ez);
    void getMaterial (int e, material_t &material);

    void calc_bmat_3D (int gp, double bmat[6][3*8]);

    void calcDistributions (void);
    void calcAverageStress (double stress_ave[6]);
    void calcAverageStrain (double strain_ave[6]);

    void writeVtu (int time_step, int elem);
    void output (int time_step, int elem, int macro_gp_global, double *MacroStrain);

};
