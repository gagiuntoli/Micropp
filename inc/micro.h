#include <vector>
#include <iostream>
#include <list>
#include "ell.h"

#define MAX_MAT_PARAM 10
#define MAX_MATS      10
#define MAX_GP_VARS   10

#define glo_elem3D(ex,ey,ez) ((ez)*(nx-1)*(ny-1) + (ey)*(nx-1) + (ex))

using namespace std;

struct MacroGp_t {
  int id;
  double *int_vars;
};

struct material_t {
  double E;
  double nu;
  double k;
  double mu;
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

    double *strain; // average strain on each element
    double *stress; // average stress on each element
    int *elem_type; // number that is changes depending on the element type
    double* int_vars; // internal variables at each Gauss point

    int NewRap_Its;
    double NewRap_Tol;

    double stress_ave[6];
    double strain_ave[6];

    ell_matrix A;
    ell_solver solver;
    double *u, *du, *b;

    Problem (int dim, int size[3], int micro_type, double *micro_params, int *mat_types, double *params);
    ~Problem (void);

    void loc_hom_Stress (int macro_id, double *MacroStrain, double *MacroStress);
    void loc_hom_Ctan (int macroGp_id, double *MacroStrain, double *MacroCtan);

    void setDisp (double *eps);

    void Assembly_A (void);
    double Assembly_b (void);

    void solve (void);
    void newtonRaphson (void);

    void getElemental_A (int ex, int ey, double (&Ae)[2*4*2*4]);
    void getElemental_A (int ex, int ey, int ez, double (&Ae)[3*8*3*8]);

    void getElemental_b (int ex, int ey, double (&be)[2*4]);
    void getElemental_b (int ex, int ey, int ez, double (&be)[3*8]);

    void getStrain (int ex, int ey, int gp, double *strain_gp);
    void getStrain (int ex, int ey, int ez, int gp, double *strain_gp);

    void getStress (int ex, int ey, int gp, double *stress_gp);
    void getStress (int ex, int ey, int ez, int gp, double *stress_gp);

    void getElemDisp (int ex, int ey, double *elem_disp);
    void getElemDisp (int ex, int ey, int ez, double *elem_disp);

    int getElemType (int ex, int ey);
    int getElemType (int ex, int ey, int ez);
    void getMaterial (int e, material_t &material);

    void calc_bmat_3D (int gp, double bmat[6][3*8]);

    void calcDistributions (void);
    void calcAverageStress (void);
    void calcAverageStrain (void);

    void writeVtu (int time_step, int elem) ;
    void output (int time_step, int elem, int macro_gp_global, double *MacroStrain);

};
