#include <vector>

class Problem {

  public:

    int npe = 4;
    int dim = 2;
    int nvoi= 3;

    std::vector<int> bc_nods, bc_y0, bc_y1, bc_x0, bc_x1;
    double xg[4][2];
    double wg[4];
    double b_mat[3][8][4];
    double dsh[4][2][4];
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

    Problem (int argc, char *argv[]);
};
