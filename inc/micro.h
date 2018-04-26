#include <vector>
#include "csr.h"
#include "ell.h"

class Problem {

  public:

    int npe = 4;
    int dim = 2;
    int nvoi= 3;

    std::vector<int> bc_nods, bc_y0, bc_y1, bc_x0, bc_x1;
    double wg[4];
    int nx, ny, nz, nn; 
    double lx, ly, lz, dx, dy, dz;
    int nelem;
    int X0Y0_nod;
    int X1Y0_nod;
    int X1Y1_nod;
    int X0Y1_nod;
    int size_tot;
    int elem_type;
    std::vector<double> strain;
    std::vector<double> stress;
    std::vector<std::vector<int> > elements;
    std::vector<double> int_vars;
    double Ef, Em;
    double Sy_f, Sy_m;

    ell_matrix A_ell;
    ell_solver solver_ell;
    double *u;
    double *du;
    double *b;

    Problem (int dim, int size[3]);
    void assembly_A (void);
    double assembly_b (void);
    void solve (void);
    void get_elemental_A (int e, double (&Ae)[8][8]);
    void get_elemental_b (int e, double (&be)[8]);
    double distance (int e);
    void write_vtu (void);

};
