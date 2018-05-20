#include <vector>
#include <list>
#include "ell.h"

struct micro_gauss_point_t {
  int micro_gp;
  double *int_vars;
};

struct macro_gauss_point {
  int macro_gp;
  std::list<micro_gauss_point_t> micro_gauss_points;
};

class Problem {

  public:

    int npe;
    int dim;
    int nvoi;

    int nx, ny, nz, nn; 
    double lx, ly, lz, dx, dy, dz;
    int nelem;
    int size_tot;

    int micro_type;

    double *strain; // average strain on each element
    double *stress; // average stress on each element
    int *elem_type; // number that is changes depending on the element type
    double* int_vars; // internal variables at each Gauss point

    double Ef, Em;
    double Sy_f, Sy_m;

    int NewRap_Its;
    double NewRap_Tol;

    double stress_ave[6];
    double strain_ave[6];

    ell_matrix A;
    ell_solver solver;
    double *u;
    double *du;
    double *b;

    bool flag_print_A;
    bool flag_print_b;
    bool flag_print_u;
    bool flag_print_du;
    bool flag_print_newton;
    bool flag_print_solver;
    bool flag_print_wrapper;

    Problem (int dim, int size[3], int cg_its, double cg_tol);
    ~Problem (void);

    void loc_hom_Stress (double *MacroStrain, double *MacroStress);
    void loc_hom_Ctan   (double *MacroStrain, double *MacroCtan);

    void setDisp (double *eps);

    void Assembly_A (void);
    double Assembly_b (void);

    void solve (void);
    void newtonRaphson (void);

    void getElemental_A (int e, double (&Ae)[3*8*3*8]);
    void getElemental_A (int ex, int ey, int ez, double (&Ae)[3*8*3*8]);

    double distance (int e);
    double distance (int ex, int ey, int ez);

    void getElemental_b (int e, double (&be)[3*8]);
    void getElemental_b (int ex, int ey, int ez, double (&be)[3*8]);

    void getStrain (int e, int gp, double *strain_gp);
    void getStrain (int ex, int ey, int ez, int gp, double *strain_gp);

    void getStress (int e, int gp, double *stress_gp);
    void getStress (int ex, int ey, int ez, int gp, double *stress_gp);

    void getElemDisp (int e, double *elem_disp);
    void getElemDisp (int ex, int ey, int ez, double *elem_disp);

    int getElemType (int e);
    int getElemType (int ex, int ey, int ez);

    void calc_bmat_3D (int gp, double bmat[6][3*8]);

    void calcDistributions (void);
    void calcAverageStress (void);
    void calcAverageStrain (void);

    void writeVtu (int time_step, int elem) ;
    void output (int time_step, int elem, int macro_gp_global, double *MacroStrain);

//  private:

};
