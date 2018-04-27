#include <vector>
#include "ell.h"

class Problem {

  public:

    int npe = 4;
    int dim = 2;
    int nvoi= 3;

    double wg[4];
    int nx, ny, nz, nn; 
    double lx, ly, lz, dx, dy, dz;
    int nelem;
    int size_tot;

    double *strain; // average strain on each element
    double *stress; // average stress on each element
    int *elem_type; // number that is changes depending on the element type
    double* int_vars; // internal variables at each Gauss point

    double Ef, Em;
    double Sy_f, Sy_m;

    ell_matrix A;
    ell_solver solver;
    double *u;
    double *du;
    double *b;

    bool flag_print_A;
    bool flag_print_b;
    bool flag_print_u;
    bool flag_print_du;

    Problem (int dim, int size[3], int cg_its, double cg_tol);
    void assembly_A (void);
    double assembly_b (void);
    void solve (void);
    void NewtonRaphson (void);
    void get_elemental_A (int e, double (&Ae)[8][8]);
    void getElemental_b (int e, double (&be)[8]);
    double distance (int e);
    void write_vtu (void);
    void getStrain (int e, int gp, double *strain_gp);
    void getStress (int e, int gp, double *stress_gp);
    void getElemDisp (int e, double *elem_disp);
    int getElemType (int e);

//  private:


};
